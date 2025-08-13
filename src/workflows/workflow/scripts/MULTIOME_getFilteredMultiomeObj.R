###############################################################################
# MULTIOME_getFilteredMultiomeObj.R
# Script to filter a multiome object, perform dimensionality reduction,
# clustering, and save results. Modular, production-ready version.
#
# Usage: Rscript MULTIOME_getFilteredMultiomeObj.R --date ... [other args]
############

# ================================
# Libraries
# ================================
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(future)
  library(ggplot2)
  library(argparse)
  library(dplyr)
  library(Seurat)
})

# parse arguments
parse_args <- function() {
  parser <- ArgumentParser(description = "Filter and process multiome object.")
  parser$add_argument("--date", type = "character", required = TRUE,
                     help = "Date for output file naming.")
  parser$add_argument("--figdir", type = "character", required = TRUE,
                     help = "Directory for saving figures (must exist).")
  parser$add_argument("--datadir", type = "character", required = TRUE,
                     help = "Directory for data files (must exist).")
  parser$add_argument("--workingdir", type = "character", required = TRUE,
                     help = "Working directory (must exist).")
  parser$add_argument("--unfilteredMultiomeFile", type = "character", required = TRUE,
                     help = "Path to unfiltered multiome .rds file.")
  parser$add_argument("--filteredMultiomeFileSaveName", type = "character", required = TRUE,
                     help = "Filename to save filtered multiome object.")
  parser$add_argument("--sampleName", type = "character", required = TRUE,
                     help = "Sample name for barcode prefix.")
  parser$add_argument("--nFeature_RNA_max", type = "double", default = 8000)
  parser$add_argument("--nFeature_RNA_min", type = "double", default = 500)
  parser$add_argument("--nFeature_ATAC_max", type = "double", default = 80000)
  parser$add_argument("--nFeature_ATAC_min", type = "double", default = 1000)
  parser$add_argument("--percent_mt_max", type = "double", default = 20)
  parser$add_argument("--TSS_enrichment_min", type = "double", default = 2)
  parser$add_argument("--outsdir", type = "character", required = FALSE, default = NULL)
  args <- parser$parse_args()
  return(args)
}

# -----------------------------
# Functions
# -----------------------------
get_metadata_filename <- function(filtered_file, date) {
  sub(".rds$", paste0("_", date, "metadata.csv"), filtered_file)
}
get_data_figure_name = function(dir_figures, figure_type, assay, date, suffix){
        figname = paste(c(figure_type, assay, date), collapse = '-')
        fullfigname = paste0(dir_figures, figname, suffix)
        return(fullfigname)  
}

#' Filter cells based on RNA/ATAC features, percent.mt, and TSS enrichment
#' @param meta_data Dataframe of cell metadata
#' @param nFeature_RNA_max, nFeature_RNA_min, nFeature_ATAC_max, nFeature_ATAC_min, percent_mt_max, TSS_enrichment_min Numeric thresholds
#' @return Filtered metadata dataframe
filter_cells <- function(meta_data, nFeature_RNA_max, nFeature_RNA_min, nFeature_ATAC_max, nFeature_ATAC_min, percent_mt_max, TSS_enrichment_min) {
  meta_data %>%
    filter(
      nFeature_ATAC < nFeature_ATAC_max & nFeature_RNA < nFeature_RNA_max &
      nFeature_ATAC > nFeature_ATAC_min & nFeature_RNA > nFeature_RNA_min &
      percent.mt < percent_mt_max & TSS.enrichment >= TSS_enrichment_min
    )
}

#' Run RNA workflow: SCTransform, PCA, UMAP, and plot depth correlation
run_rna_workflow <- function(sample, figures_dir, sample_name, date) {
  DefaultAssay(sample) <- "RNA"
  sample <- SCTransform(sample, verbose = TRUE)
  sample <- RunPCA(sample)
  plot_name <- get_data_figure_name(figures_dir, paste0("DepthCor-pca-", sample_name), "RNA", date, ".pdf")
  pdf(plot_name)
  DepthCor(sample, assay = 'RNA', reduction = "pca")
  dev.off()
  sample <- RunUMAP(sample, dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
  cat("Completed RNA UMAP.\n")
  return(sample)
}

#' Run ATAC workflow: TFIDF, SVD, UMAP, and plot depth correlation
run_atac_workflow <- function(sample, figures_dir, sample_name, date) {
  DefaultAssay(sample) <- "ATAC"
  sample <- RunTFIDF(sample)
  sample <- FindTopFeatures(sample, min.cutoff = 'q0')
  sample <- RunSVD(sample)
  depthcors <- DepthCorLSI(sample, assay = 'ATAC', reduction = "lsi")
  dims <- which(!(abs(depthcors$counts) > 0.7))
  cat(sprintf("ATAC depth correlation: %s\n", paste(depthcors$counts, collapse = ", ")))
  cat(sprintf("Keeping LSI dims: %s\n", paste(dims, collapse = ", ")))
  plot_name <- get_data_figure_name(figures_dir, paste0("DepthCor-lsi-", sample_name), "ATAC", date, ".pdf")
  pdf(plot_name)
  DepthCor(sample, assay = 'ATAC', reduction = "lsi")
  dev.off()
  sample <- RunUMAP(sample, reduction = 'lsi', dims = dims, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
  cat("Completed ATAC UMAP.\n")
  return(list(sample = sample, dims = dims))
}

#' Run multimodal neighbors, clustering, and plot UMAPs
run_multimodal_and_plot <- function(sample, figures_dir, sample_name, date, dims) {
  cat("Running multimodal neighbors and clustering...\n")
  sample <- FindMultiModalNeighbors(sample, reduction.list = list("pca", "lsi"), dims.list = list(1:50, dims))
  sample <- RunUMAP(sample, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
  sample <- FindClusters(sample, graph.name = "wsnn", algorithm = 3, verbose = TRUE, resolution = c(1.0, 0.4, 0.8, 0.6))
  plot_name <- get_data_figure_name(figures_dir, paste0("UMAPs-", sample_name), "MULTIOME", date, ".png")
  p1 <- DimPlot(sample, reduction = "umap.rna", group.by = "wsnn_res.0.8", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle(paste0(sample_name, " RNA"))
  p2 <- DimPlot(sample, reduction = "umap.atac", group.by = "wsnn_res.0.8", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle(paste0(sample_name, " ATAC"))
  p3 <- DimPlot(sample, reduction = "wnn.umap", group.by = "wsnn_res.0.8", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle(paste0(sample_name, " WNN"))
  p4 <- DimPlot(sample, reduction = "umap.rna", group.by = "wsnn_res.0.6", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle(paste0(sample_name, " RNA"))
  p5 <- DimPlot(sample, reduction = "umap.atac", group.by = "wsnn_res.0.6", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle(paste0(sample_name, " ATAC"))
  p6 <- DimPlot(sample, reduction = "wnn.umap", group.by = "wsnn_res.0.6", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle(paste0(sample_name, " WNN"))
  p7 <- DimPlot(sample, reduction = "umap.rna", group.by = "wsnn_res.0.4", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle(paste0(sample_name, " RNA"))
  p8 <- DimPlot(sample, reduction = "umap.atac", group.by = "wsnn_res.0.4", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle(paste0(sample_name, " ATAC"))
  p9 <- DimPlot(sample, reduction = "wnn.umap", group.by = "wsnn_res.0.4", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle(paste0(sample_name, " WNN"))
  # Combine and save plots
  p_all <- p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
  ggsave(plot_name, plot = p_all, width = 15, height = 15)
  cat("Saved UMAP plots.\n")
  return(sample)
}


# main fucnction to run per sample

main <- function() {
  args <- parse_args()
  setwd(args$workingdir)

  # Read sample
  if (!file.exists(args$unfilteredMultiomeFile)) {
    cat(sprintf("Unfiltered multiome file not found: %s\n", args$unfilteredMultiomeFile))
    stop("Input file missing.")
  }
  sample <- readRDS(args$unfilteredMultiomeFile)
  cat(sprintf("Loaded unfiltered sample. Dimensions: %s\n", paste(dim(sample), collapse = ", ")))

  # Filter cells
  filtered_meta <- filter_cells(
    sample@meta.data,
    args$nFeature_RNA_max, args$nFeature_RNA_min,
    args$nFeature_ATAC_max, args$nFeature_ATAC_min,
    args$percent_mt_max, args$TSS_enrichment_min
  )
  sample <- subset(sample, cells = rownames(filtered_meta))
  cat(sprintf("Filtered sample. Dimensions: %s\n", paste(dim(sample), collapse = ", ")))

  # RNA workflow
  sample <- run_rna_workflow(sample, args$figdir, args$sampleName, args$date)

  # ATAC workflow
  atac_result <- run_atac_workflow(sample, args$figdir, args$sampleName, args$date)
  sample <- atac_result$sample
  dims <- atac_result$dims

  # Multimodal neighbors, clustering, and UMAP plots
  sample <- run_multimodal_and_plot(sample, args$figdir, args$sampleName, args$date, dims)

  # Save results
  metadata_filename <- get_metadata_filename(args$filteredMultiomeFileSaveName, args$date)
  cat(sprintf("Saving filtered object to: %s\n", args$filteredMultiomeFileSaveName))
  saveRDS(sample, args$filteredMultiomeFileSaveName)
  write.csv(sample@meta.data, metadata_filename)
  cat(sprintf("Saved metadata to: %s\n", metadata_filename))
  cat("Processing complete.\n")
}

# Run it inside sys nframe so we can use the rest of the script without having to parse args 
if (sys.nframe() == 0) {
  tryCatch(
    main(),
    error = function(e) {
      cat(e$message, "\n")
      quit(status = 1)
    }
  )
}
