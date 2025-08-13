###############################################################################
# MULTIOME_getUnfilteredMultiomeObj.R
# Author: Tina Keshavarzian
# Date: 2021-07-21 
# Description: Script to merge single-sample RNA and ATAC data into an unfiltered multiome Seurat object.
############

# ================================
# Libraries
# ================================
suppressPackageStartupMessages({
  library(Signac)
  library(Seurat)
  library(GenomicRanges)
  library(future)
  library(ggplot2)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(EnsDb.Hsapiens.v86)
  library(argparse)
})


# ================================
# Argument Parsing
# ================================
parse_args <- function() {
  parser <- ArgumentParser(description = 'Merge RNA and ATAC data into an unfiltered multiome Seurat object.')
  parser$add_argument('--date', type = 'character', required = TRUE, help = 'Date for output filenames.')
  parser$add_argument('--figdir', type = 'character', required = TRUE, help = 'Directory for figures (must exist).')
  parser$add_argument('--datadir', type = 'character', required = TRUE, help = 'Directory for data (must exist).')
  parser$add_argument('--workingdir', type = 'character', required = TRUE, help = 'Working directory (must exist).')
  parser$add_argument('--chromatinObjectFile', type = 'character', required = TRUE, help = 'RDS file for ATAC assay.')
  parser$add_argument('--h5File', type = 'character', required = TRUE, help = '10x h5 file for RNA assay.')
  parser$add_argument('--sampleName', type = 'character', required = TRUE, help = 'Sample name for barcode prefix.')
  args <- parser$parse_args()
  return(args)
}

# ================================
# Utility Functions
# ================================
validate_file_exists <- function(filepath, description) {
  if (!file.exists(filepath)) {
    cat(paste('File not found:', filepath, '(', description, ')'))
    stop(paste('File not found:', filepath, '(', description, ')'))
  }
}

validate_dir_exists <- function(dirpath, description) {
  if (!dir.exists(dirpath)) {
    cat(paste('Directory not found:', dirpath, '(', description, ')'))
    stop(paste('Directory not found:', dirpath, '(', description, ')'))
  }
}

save_seurat_object <- function(obj, rds_path, meta_csv_path) {
  cat(paste('Saving Seurat object to', rds_path))
  saveRDS(obj, rds_path)
  cat(paste('Saving metadata to', meta_csv_path))
  write.csv(obj@meta.data, meta_csv_path)
}

create_multiome_object <- function(rna_h5, atac_rds, sample_name) {
  # Read RNA data
  cat('Reading RNA h5 file...')
  h5_sample <- Read10X_h5(rna_h5)
  rna_sample <- h5_sample$`Gene Expression`
  rna_seurat <- CreateSeuratObject(counts = rna_sample)
  rna_seurat <- RenameCells(rna_seurat, add.cell.id = sample_name)
  cat('RNA Seurat object created.')

  # Read ATAC data
  cat('Reading ATAC RDS file...')
  atac_seurat <- readRDS(atac_rds)
  atac_seurat <- RenameCells(atac_seurat, add.cell.id = sample_name)
  cat('ATAC Seurat object loaded.')

  # Merge into multiome object
  cat('Merging RNA and ATAC into multiome object...')
  multiome <- rna_seurat
  multiome[['ATAC']] <- atac_seurat[['ATAC']]
  multiome[['percent.mt']] <- PercentageFeatureSet(multiome, pattern = '^MT-')
  multiome@meta.data <- cbind(multiome@meta.data, atac_seurat@meta.data)
  if (all(c('atac_peak_region_fragments', 'atac_fragments') %in% colnames(multiome@meta.data))) {
    multiome$frip <- multiome@meta.data$atac_peak_region_fragments / multiome@meta.data$atac_fragments
  } else {
    cat('FRIP calculation skipped: required columns missing in metadata.')
    multiome$frip <- NA
  }
  multiome$name <- rep(sample_name, ncol(multiome))
  cat('Multiome object created.')
  return(multiome)
}

# ================================
# Main
# ================================
main <- function() {
  args <- parse_args()

  # Validate directories and files
  validate_dir_exists(args$figdir, 'Figures directory')
  validate_dir_exists(args$datadir, 'Data directory')
  validate_dir_exists(args$workingdir, 'Working directory')
  validate_file_exists(args$chromatinObjectFile, 'ATAC RDS file')
  validate_file_exists(args$h5File, 'RNA h5 file')

  setwd(args$workingdir)

  # Output filenames
  multiome_rds <- file.path(args$datadir, paste0('multiome_object_clustered_', args$sampleName, '_', args$date, '.rds'))
  multiome_unfiltered_rds <- file.path(args$datadir, paste0('multiome_object_unfiltered_', args$sampleName, '_', args$date, '.rds'))
  multiome_unfiltered_meta_csv <- file.path(args$datadir, paste0('multiome_object_unfiltered_', args$sampleName, '_', args$date, 'metadata.csv'))

  # Create multiome object
  multiome_obj <- create_multiome_object(args$h5File, args$chromatinObjectFile, args$sampleName)

  # Save outputs
  save_seurat_object(multiome_obj, multiome_unfiltered_rds, multiome_unfiltered_meta_csv)
  cat('All done!')
}

if (sys.nframe() == 0) {
  tryCatch(
    main(),
    error = function(e) {
      cat(paste('Script failed:', e$message))
      stop(e)
    }
  )
}
