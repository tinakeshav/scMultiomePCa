# merging multiome objects

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(future)
  library(ggplot2)
  library(argparse)
})

# ================================
# Utility Functions
# ================================

#' Log a message with timestamp; merged takes a while so I'm time stamping it
#' @param ... Message(s) to print
log_message <- function(...) {
  cat(sprintf("[%s] ", Sys.time()), ... , "\n")
}

#' Check if required directories exist
#' @param dirs Character vector of directory paths
check_directories_exist <- function(dirs) {
  for (dir in dirs) {
    if (!dir.exists(dir)) {
      stop(sprintf("Directory does not exist: %s", dir))
    }
  }
}

#' Read multiple RDS files into a list
#' @param file_paths Character vector of file paths
#' @return List of loaded objects
read_multiome_files <- function(file_paths) {
  lapply(file_paths, function(f) {
    if (!file.exists(f)) stop(sprintf("File not found: %s", f))
    readRDS(f)
  })
}

#' Save merged object and metadata
#' @param object The merged Seurat object
#' @param rds_path Path to save RDS
#' @param metadata_path Path to save metadata CSV
save_merged_outputs <- function(object, rds_path, metadata_path) {
  saveRDS(object, rds_path)
  write.csv(object@meta.data, metadata_path)
}

# ================================
# Main Analysis Pipeline
# ================================

main <- function() {
  # Argument parsing
  parser <- ArgumentParser(description = "Merge Seurat multiome objects and process.")
  parser$add_argument("--date", type = "character", required = TRUE,
                      help = "Date script is running; used in filenames.")
  parser$add_argument("--figdir", type = "character", required = TRUE,
                      help = "Figures directory; must already exist.")
  parser$add_argument("--datadir", type = "character", required = TRUE,
                      help = "Data directory; must already exist.")
  parser$add_argument("--workingdir", type = "character", required = TRUE,
                      help = "Working directory; must already exist.")
  parser$add_argument("--multiomeFiles", type = "character", nargs = '+', required = TRUE,
                      help = "Filtered multiome files to read.")
  parser$add_argument("--mergedMultiomeFileName", type = "character", required = TRUE,
                      help = "Merged multiome filename to save as.")

  args <- parser$parse_args()
  log_message("Parsed arguments:")
  print(args)

  # Set working directory and check directories
  check_directories_exist(c(args$figdir, args$datadir, args$workingdir))
  setwd(args$workingdir)
  log_message("Set working directory to", args$workingdir)

  # Source utility functions
  source('scripts/sc_utilFunctions.R')

  # Prepare output filenames
  merged_metadata_name <- sub('.rds$', paste0('_', args$date, 'metadata.csv'), args$mergedMultiomeFileName)
  name_to_keep <- 'merged'

  # Read and merge objects
  log_message("Reading multiome files...")
  multiome_objects <- read_multiome_files(args$multiomeFiles)
  log_message("Read", length(multiome_objects), "objects.")

  log_message("Merging Seurat objects...")
  sample <- merge_seurat_list(multiome_objects)
  log_message("Merged all objects.")

  log_message("Sample dimensions:", paste(dim(sample), collapse = ' x '))

  # RNA analysis
  log_message("Running RNA analysis...")
  DefaultAssay(sample) <- "RNA"
  sample <- SCTransform(sample, verbose = TRUE)
  sample <- RunPCA(sample)

  plot_depthcor_pca_name <- get_data_figure_name(args$figdir, paste0('DepthCor-pca-', name_to_keep), 'RNA', args$date, '.pdf')
  pdf(plot_depthcor_pca_name)
  DepthCor(sample, assay = 'RNA', reduction = "pca")
  dev.off()

  sample <- RunUMAP(sample, dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
  log_message("Ran UMAP on RNA.")

  # ATAC analysis
  log_message("Running ATAC analysis...")
  DefaultAssay(sample) <- "ATAC"
  sample <- RunTFIDF(sample)
  sample <- FindTopFeatures(sample, min.cutoff = 'q0')
  sample <- RunSVD(sample)
  depthcors <- DepthCorLSI(sample, assay = 'ATAC', reduction = "lsi")
  dims <- which(!(abs(depthcors$counts) > 0.7))
  log_message("DepthCorLSI dimensions kept:", paste(dims, collapse = ','))

  plot_depthcor_svd_name <- get_data_figure_name(args$figdir, paste0('DepthCor-lsi-', name_to_keep), 'ATAC', args$date, '.pdf')
  pdf(plot_depthcor_svd_name)
  DepthCor(sample, assay = 'ATAC', reduction = "lsi")
  dev.off()

  sample <- RunUMAP(sample, reduction = 'lsi', dims = dims, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
  log_message("Ran UMAP on ATAC.")

  # Multi-modal neighbors and clustering
  log_message("Running multi-modal neighbors and clustering...")
  sample <- FindMultiModalNeighbors(sample, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
  sample <- RunUMAP(sample, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
  sample <- FindClusters(sample, graph.name = "wsnn", algorithm = 3, verbose = TRUE, resolution = c(1.0, 0.4, 0.8, 0.6))
  log_message("Clustering complete. Meta data columns:", paste(colnames(sample@meta.data), collapse = ', '))

  # Plot UMAPs
  log_message("Plotting UMAPs...")
  plot_umaps_name <- get_data_figure_name(args$figdir, paste0('UMAPs-', name_to_keep), 'MULTIOME', args$date, '.png')
  p1 <- DimPlot(sample, reduction = "umap.rna", group.by = "wsnn_res.0.8", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle(paste0(name_to_keep, " RNA"))
  p2 <- DimPlot(sample, reduction = "umap.atac", group.by = "wsnn_res.0.8", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle(paste0(name_to_keep, " ATAC"))
  p3 <- DimPlot(sample, reduction = "wnn.umap", group.by = "wsnn_res.0.8", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle(paste0(name_to_keep, " WNN"))
  p4 <- DimPlot(sample, reduction = "umap.rna", group.by = "wsnn_res.0.6", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle(paste0(name_to_keep, " RNA"))
  p5 <- DimPlot(sample, reduction = "umap.atac", group.by = "wsnn_res.0.6", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle(paste0(name_to_keep, " ATAC"))
  p6 <- DimPlot(sample, reduction = "wnn.umap", group.by = "wsnn_res.0.6", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle(paste0(name_to_keep, " WNN"))
  p7 <- DimPlot(sample, reduction = "umap.rna", group.by = "wsnn_res.0.4", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle(paste0(name_to_keep, " RNA"))
  p8 <- DimPlot(sample, reduction = "umap.atac", group.by = "wsnn_res.0.4", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle(paste0(name_to_keep, " ATAC"))
  p9 <- DimPlot(sample, reduction = "wnn.umap", group.by = "wsnn_res.0.4", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle(paste0(name_to_keep, " WNN"))
  p_all <- p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
  ggsave(plot_umaps_name, p_all, width = 15, height = 15)

  plot_umaps_byname <- get_data_figure_name(args$figdir, paste0('UMAPs-byName', name_to_keep), 'MULTIOME', args$date, '.png')
  p1n <- DimPlot(sample, reduction = "umap.rna", group.by = "name", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle(paste0(name_to_keep, " RNA"))
  p2n <- DimPlot(sample, reduction = "umap.atac", group.by = "name", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle(paste0(name_to_keep, " ATAC"))
  p3n <- DimPlot(sample, reduction = "wnn.umap", group.by = "name", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle(paste0(name_to_keep, " WNN"))
  p_byname <- p1n + p2n + p3n & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
  ggsave(plot_umaps_byname, p_byname, width = 15, height = 5)

  # Save outputs
  log_message("Saving merged object and metadata...")
  save_merged_outputs(sample, args$mergedMultiomeFileName, merged_metadata_name)
  log_message("Saved merged object to", args$mergedMultiomeFileName)
  log_message("Saved metadata to", merged_metadata_name)
}

# ================================
# Run main if script is executed 
# otherwise just sources the script 
# ================================

if (sys.nframe() == 0) {
  tryCatch(
    main(),
    error = function(e) {
      cat(sprintf("[ERROR %s] %s\n", Sys.time(), e$message))
      quit(status = 1)
    }
  )
}
