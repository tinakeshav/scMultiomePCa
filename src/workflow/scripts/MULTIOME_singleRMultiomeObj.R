###############################################################################
# MULTIOME_singleRMultiomeObj.R
# Description: Run SingleR with multiple reference atlases
############

# ================================
# Libraries
# ================================
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(future)
  library(ggplot2)
  library(argparse)
  library(Signac)
  library(Seurat)
  library(SingleR)
  library(celldex)
  library(data.table)
})

# ================================
# Argument Parsing
# ================================
parser <- ArgumentParser(description = 'Annotate multiome Seurat object using SingleR and reference atlases.')
parser$add_argument('--date', type = 'character', required = TRUE,
                    help = 'Date script is running; used for output filenames.')
parser$add_argument('--figdir', type = 'character', required = TRUE,
                    help = 'Directory for figures (must exist).')
parser$add_argument('--datadir', type = 'character', required = TRUE,
                    help = 'Directory for data (must exist).')
parser$add_argument('--workingdir', type = 'character', required = TRUE,
                    help = 'Working directory (must exist).')
parser$add_argument('--filteredMultiomeFile', type = 'character', required = TRUE,
                    help = 'Path to filtered multiome RDS file.')
parser$add_argument('--sampleName', type = 'character', required = TRUE,
                    help = 'Sample name, added to barcode prefix.')

args <- parser$parse_args()
setwd(args$workingdir)
cat('Arguments parsed and working directory set.\n')

# ================================
# Constants and Parameters
# ================================
DATE <- args$date
DATA_DIR <- args$datadir
FIGURES_DIR <- args$figdir
MULTIOME_OBJ_PATH <- args$filteredMultiomeFile
SAMPLE_NAME <- args$sampleName
CLUSTERS_TO_CHECK <- c('wsnn_res.0.4', 'wsnn_res.0.6', 'wsnn_res.0.8', 'wsnn_res.1')
REDUCTION_NAMES_TO_PLOT <- c('umap.rna', 'umap.atac', 'wnn.umap')

# ================================
# Helper Functions
# ================================

#' Log-normalize a SingleCellExperiment atlas and set rownames if needed
#' @param scrna_atlas SingleCellExperiment object
#' @param label Column name for cell type labels
#' @param rowdata_symbol Logical, whether to set rownames from rowData
#' @param symbol Column name in rowData for gene symbols
#' @return Log-normalized SingleCellExperiment
scrna_analysis <- function(scrna_atlas, label, rowdata_symbol = FALSE, symbol = 'symbol') {
  scrna_atlas <- scrna_atlas[, !is.na(scrna_atlas[[label]])]
  scrna_atlas <- scrna_atlas[, colSums(assays(scrna_atlas)$counts) > 0] # Remove libraries with no counts
  scrna_atlas <- scuttle::logNormCounts(scrna_atlas)
  if (rowdata_symbol) {
    cat('Setting rownames from rowData symbol column.\n')
    rownames(scrna_atlas) <- rowData(scrna_atlas)[, symbol]
  }
  return(scrna_atlas)
}

#' Run SingleR annotation and add metadata to Seurat object
#' @param seurat_obj Seurat object
#' @param refs Named list of reference atlases
#' @param ref_labels List of label vectors for each reference
#' @param clusters_to_check Character vector of cluster column names
#' @param reduction_names_to_plot Character vector of reduction names
#' @param sample_name Sample name string
#' @param date Date string
#' @return Annotated Seurat object
run_singleR_and_annotate <- function(seurat_obj, refs, ref_labels, clusters_to_check,
                                     reduction_names_to_plot, sample_name, date, data_dir, figures_dir) {
  cat('Running SingleR annotation with combined references.\n')
  pred <- SingleR(GetAssayData(seurat_obj, assay = 'RNA', slot = 'data'),
                  ref = refs, labels = ref_labels)
  save_path <- file.path(data_dir, paste0('multiome_object_filtered_singleR_combined_int_pred_',
                                          sample_name, '_clusters_', date, '.rds'))
  saveRDS(pred, save_path)
  cat(sprintf('SingleR predictions saved to %s\n', save_path))

  singleR_datasets <- names(refs)
  for (dataset in singleR_datasets) {
    metadata_name <- paste0('bulk', dataset, '_prunedlabels')
    cat(sprintf('Adding metadata: %s\n', metadata_name))
    seurat_obj <- AddMetaData(seurat_obj,
                             pred[['orig.results']][[dataset]]$pruned.labels,
                             col.name = metadata_name)
    # Plotting
    plot_singleR(seurat_sample = seurat_obj,
                 singleR_pred = pred,
                 singleR_pred_dataset_name = dataset,
                 clusters_to_check = clusters_to_check,
                 name_to_keep = sample_name,
                 date = date)
    for (reduction_name in reduction_names_to_plot) {
      plot_singleR_umap(seurat_sample = seurat_obj,
                        metadata_colname = metadata_name,
                        reduction_name_to_plot = reduction_name,
                        name_to_keep = sample_name,
                        date = date)
    }
  }
  return(seurat_obj)
}

#' Run SingleR annotation for scBRCA reference and add metadata
#' @param seurat_obj Seurat object
#' @param sc_brca_bach scBRCA reference object
#' @param clusters_to_check Character vector of cluster column names
#' @param reduction_names_to_plot Character vector of reduction names
#' @param sample_name Sample name string
#' @param date Date string
#' @param data_dir Data directory
#' @return Annotated Seurat object
run_singleR_scbrca <- function(seurat_obj, sc_brca_bach, clusters_to_check,
                              reduction_names_to_plot, sample_name, date, data_dir) {
  cat('Running SingleR annotation with scBRCA reference.\n')
  pred_sc_brca <- SingleR(GetAssayData(seurat_obj, assay = 'RNA', slot = 'data'),
                         ref = list(scBRCA = sc_brca_bach),
                         labels = list(sc_brca_bach$label.main))
  save_path <- file.path(data_dir, paste0('multiome_object_filtered_singleR_combined_int_pred_sc_BRCA_',
                                          sample_name, '_clusters_', date, '.rds'))
  saveRDS(pred_sc_brca, save_path)
  cat(sprintf('scBRCA SingleR predictions saved to %s\n', save_path))

  metadata_name <- 'scbrca_prunedlabels'
  seurat_obj <- AddMetaData(seurat_obj,
                           pred_sc_brca[['orig.results']][['scBRCA']]$pruned.labels,
                           col.name = metadata_name)
  # Plotting
  plot_singleR(seurat_sample = seurat_obj,
               singleR_pred = pred_sc_brca,
               singleR_pred_dataset_name = 'scBRCA',
               clusters_to_check = clusters_to_check,
               name_to_keep = sample_name,
               date = date)
  for (reduction_name in reduction_names_to_plot) {
    plot_singleR_umap(seurat_sample = seurat_obj,
                      metadata_colname = metadata_name,
                      reduction_name_to_plot = reduction_name,
                      name_to_keep = sample_name,
                      date = date)
  }
  return(seurat_obj)
}

# ================================
# Main Workflow
# ================================

tryCatch({
  cat(sprintf('Loading multiome Seurat object from %s\n', MULTIOME_OBJ_PATH))
  sample <- readRDS(MULTIOME_OBJ_PATH)

  # Load reference atlases
  cat('Loading reference atlases.\n')
  hpca <- readRDS('~/common/HumanPrimaryCellAtlasData.rds')
  bpe <- readRDS('~/common/BlueprintEncodeData.rds')
  imm <- readRDS('~/common/DatabaseImmuneCellExpressionData.rds')

  # Prepare reference atlases
  hpca$label.main <- paste0('HPCA.', hpca$label.main)
  bpe$label.main <- paste0('BPE.', bpe$label.main)
  imm$label.main <- paste0('IMM.', imm$label.main)

  refs <- list(BPE = bpe, HPCA = hpca, IMM = imm)
  ref_labels <- list(bpe$label.main, hpca$label.main, imm$label.main)

  # Annotate with SingleR using combined references
  sample <- run_singleR_and_annotate(
    seurat_obj = sample,
    refs = refs,
    ref_labels = ref_labels,
    clusters_to_check = CLUSTERS_TO_CHECK,
    reduction_names_to_plot = REDUCTION_NAMES_TO_PLOT,
    sample_name = SAMPLE_NAME,
    date = DATE,
    data_dir = DATA_DIR,
    figures_dir = FIGURES_DIR
  )

  # Load and prepare scBRCA reference
  cat('Loading and preparing scBRCA reference.\n')
  sc_brca_bach <- readRDS('/cluster/projects/lupiengroup/People/tina/bca/data/processed/BRCA1_SCE.rds')
  sc_brca_bach$label.main <- sc_brca_bach$Groups
  rownames(sc_brca_bach) <- toupper(rownames(sc_brca_bach))
  sc_brca_bach <- scrna_analysis(sc_brca_bach, 'Groups', rowdata_symbol = FALSE)

  # Annotate with SingleR using scBRCA
  sample <- run_singleR_scbrca(
    seurat_obj = sample,
    sc_brca_bach = sc_brca_bach,
    clusters_to_check = CLUSTERS_TO_CHECK,
    reduction_names_to_plot = REDUCTION_NAMES_TO_PLOT,
    sample_name = SAMPLE_NAME,
    date = DATE,
    data_dir = DATA_DIR
  )

  # Save metadata
  metadata_csv_path <- sub('.rds$', paste0('_', DATE, 'singleR_metadata.csv'), MULTIOME_OBJ_PATH)
  write.csv(sample@meta.data, metadata_csv_path)
  cat(sprintf('Metadata saved to %s\n', metadata_csv_path))

}, error = function(e) {
  cat(sprintf('Error occurred: %s\n', e$message))
  stop(e)
})
