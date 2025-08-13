###############################################################################
# MULTIOME_ATAC_featureMatrix_merged_peaks.R
# Script to call intra-patient peaks and compute fragment counts
###############################################################################

# -----------------------------
# Libraries
# -----------------------------
suppressPackageStartupMessages({
  library(argparse)
  library(Seurat)
  library(Signac)
  library(EnsDb.Hsapiens.v86)
  library(GenomicRanges)
})

# -----------------------------
# Functions
# -----------------------------

#' Log a message with timestamp (uses cat)
#' @param msg Message to log
log_message <- function(msg) {
  cat(sprintf("[%s] %s\n", Sys.time(), msg))
}

#' Read and validate a file
#' @param file_path Path to file
#' @param read_fun Function to read file
#' @param ... Additional arguments to read_fun
#' @return Loaded object
read_and_validate <- function(file_path, read_fun, ...) {
  if (!file.exists(file_path)) {
    stop(sprintf("File not found: %s", file_path))
  }
  log_message(sprintf("Reading file: %s", file_path))
  read_fun(file_path, ...)
}

#' Save an RDS object with logging
#' @param obj Object to save
#' @param file_path Path to save
save_rds_with_log <- function(obj, file_path) {
  log_message(sprintf("Saving RDS: %s", file_path))
  saveRDS(obj, file_path)
}

#' Process a single sample: create fragment and chromatin assay objects
#' @param sample_name Name of the sample
#' @param seurat_obj Seurat object
#' @param frag_paths List of fragment paths
#' @param peaks Merged peaks object
#' @param output_fragment_prefix Output file prefix
#' @return List with fragment matrix and chromatin assay
process_sample <- function(sample_name, seurat_obj, frag_paths, peaks, output_fragment_prefix) {
  log_message(sprintf("Processing sample: %s", sample_name))
  sample_obj <- subset(seurat_obj, name == sample_name)
  frag_path <- frag_paths[grepl(paste0(sample_name, '/'), frag_paths)][[1]]
  log_message(sprintf("Fragment path: %s", frag_path))

  fragments <- CreateFragmentObject(path = frag_path,  
                                    cells = sample_obj$gex_barcode)

  fragment_matrix <- FeatureMatrix(
    fragments,
    cells = names(fragments@cells),
    features = peaks,
    verbose = TRUE
  )

  frag_rds_path <- paste0(output_fragment_prefix, '_', sample_name, '_Divided.rds')
  save_rds_with_log(fragment_matrix, frag_rds_path)

  chromatin_assay <- CreateChromatinAssay(fragment_matrix, fragments = fragments)
  chromatin_rds_path <- paste0(output_fragment_prefix, '_', sample_name, '_Divided.rds')
  save_rds_with_log(chromatin_assay, chromatin_rds_path)

  list(fragment_matrix = fragment_matrix, chromatin_assay = chromatin_assay)
}

#' Main workflow function
main <- function() {
  # -----------------------------
  # Argument Parsing
  # -----------------------------
  parser <- ArgumentParser()
  parser$add_argument("--multiomeobjname", type = 'character', required = TRUE, help = "Multiome object to load.")
  parser$add_argument("--multiomeobjpeaksRDS", type = 'character', required = TRUE, help = "Merged peaks, RDS format.")
  parser$add_argument("--metadatatoattach", type = 'character', required = TRUE, help = "Metadata file to attach.")
  parser$add_argument("--date", type = 'character', required = TRUE, help = "Date script is running; used in saving object names.")
  parser$add_argument("--output_fragment_fname", type = 'character', required = TRUE, help = "Output file name prefix for fragments.")
  parser$add_argument("--output_multiome_fname", type = 'character', required = TRUE, help = "Output file name for new multiome object.")
  parser$add_argument("--output_txtfile", type = 'character', required = TRUE, help = "Output txt file for completion stamp.")
  args <- parser$parse_args()
  print(args)

  # -----------------------------
  # Load Data
  # -----------------------------
  metadata <- read_and_validate(args$metadatatoattach, read.csv, header = TRUE, row.names = 1)
  seurat_obj <- read_and_validate(args$multiomeobjname, readRDS)
  seurat_obj <- AddMetaData(seurat_obj, metadata)
  DefaultAssay(seurat_obj) <- 'ATACv3'
  peaks <- read_and_validate(args$multiomeobjpeaksRDS, readRDS)

  # Obtain fragment paths
  frag_paths <- lapply(Fragments(seurat_obj), function(x) x@path)

  # -----------------------------
  # Process Each Sample
  # -----------------------------
  frag_obj_list <- list()
  chromatin_obj_list <- list()
  for (sample_name in unique(seurat_obj@meta.data$name)) {
    res <- tryCatch({
      process_sample(sample_name, seurat_obj, frag_paths, peaks, args$output_fragment_fname)
    }, error = function(e) {
      log_message(sprintf("Error processing sample %s: %s", sample_name, e$message))
      NULL
    })
    if (!is.null(res)) {
      frag_obj_list[[sample_name]] <- res$fragment_matrix
      chromatin_obj_list[[sample_name]] <- res$chromatin_assay
    }
  }

  # Save lists
  save_rds_with_log(frag_obj_list, paste0(args$output_fragment_fname, '_DividedNamedList.rds'))
  save_rds_with_log(chromatin_obj_list, paste0(args$output_fragment_fname, '_DividedNamedList.rds'))

  # -----------------------------
  # Merge Chromatin Assays and Save Multiome Object
  # -----------------------------
  if (length(chromatin_obj_list) == 0) {
    stop("No chromatin assays created. Exiting.")
  }
  chromatin_assay_merged <- do.call(merge, chromatin_obj_list)
  log_message("Merged chromatin assays.")
  seurat_obj[['intrapat_peaks']] <- chromatin_assay_merged
  save_rds_with_log(seurat_obj, args$output_multiome_fname)
  save_rds_with_log(chromatin_assay_merged, args$output_fragment_fname)

  # -----------------------------
  # Completion Stamp
  # -----------------------------
  log_message('Stamp of success, so snakemake knows!')
  writeLines(args$output_fragment_fname, args$output_txtfile)
}

# -----------------------------
# Run main if script is executed
# -----------------------------
if (sys.nframe() == 0) {
  tryCatch({
    main()
  }, error = function(e) {
    log_message(sprintf("Script failed: %s", e$message))
    quit(status = 1)
  })
}