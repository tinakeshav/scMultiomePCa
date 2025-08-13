###############################################################################
# ATAC_getFragmentObj_inputargs_binnedgenome.R
#   Script to generate a Seurat Chromatin Object from binned genome intervals
#   and fragment files, with optional TSS enrichment and nucleosome signal scoring.
#
# Usage:
#   Rscript ATAC_getFragmentObj_inputargs_binnedgenome.R --args ...
# Author: Tina Keshavarzian
###############################################################################

# ================================
# Libraries
# ================================
suppressPackageStartupMessages({
  library(dplyr)
  library(Signac)
  library(Seurat)
  library(GenomicRanges)
  library(future)
  library(ggplot2)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(EnsDb.Hsapiens.v86)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(argparse)
})

# ================================
# Utility Functions
# ================================

#' Read binned genome intervals from .bed or .rds file
#' @param filepath Path to .bed or .rds file
#' @return GRanges object
read_binned_genome <- function(filepath) {
  if (endsWith(filepath, ".bed")) {
    cat(sprintf("Reading binned genome from BED file: %s\n", filepath))
    df <- tryCatch({
      read.table(filepath, col.names = c("chr", "start", "end"))
    }, error = function(e) {
      cat(sprintf("Failed to read BED file: %s\n", e$message))
      stop(e)
    })
    return(makeGRangesFromDataFrame(df))
  } else if (endsWith(filepath, ".rds")) {
    cat(sprintf("Reading binned genome from RDS file: %s\n", filepath))
    return(readRDS(filepath))
  } else {
    cat(sprintf("Unsupported file type for binned genome: %s\n", filepath))
    stop("Unsupported file type for binned genome.")
  }
}

#' Read and filter barcodes file
#' @param barcodes_file Path to barcodes CSV
#' @return Data frame of filtered barcodes
read_barcodes <- function(barcodes_file) {
  cat(sprintf("Reading barcodes from: %s\n", barcodes_file))
  barcodes <- tryCatch({
    read.csv(barcodes_file)
  }, error = function(e) {
    cat(sprintf("Failed to read barcodes file: %s\n", e$message))
    stop(e)
  })
  barcodes <- barcodes[barcodes$is_cell == 1, ]
  rownames(barcodes) <- barcodes$barcode
  cat(sprintf("Barcodes dimension after filtering: %s\n", paste(dim(barcodes), collapse = ", ")))
  return(barcodes)
}

#' Create Seurat Chromatin Object
#' @param fragment_file Path to fragment file
#' @param barcodes Data frame of barcodes
#' @param binned_genome GRanges object
#' @return Seurat object
create_chromatin_object <- function(fragment_file, barcodes, binned_genome) {
  cat("Creating FragmentObject...\n")
  fragobj <- CreateFragmentObject(path = fragment_file, cells = barcodes$barcode)
  cat("Counting fragments over intervals...\n")
  counts <- FeatureMatrix(fragments = fragobj, features = binned_genome, cells = barcodes$barcode)
  cat("Creating ChromatinAssay...\n")
  atac_assay <- CreateChromatinAssay(counts = counts, fragments = fragobj)
  cat("Creating Seurat object...\n")
  chromatin_obj <- CreateSeuratObject(counts = atac_assay, assay = "ATAC", meta.data = barcodes)
  return(chromatin_obj)
}

#' Annotate and score Seurat Chromatin Object
#' @param chromatin_obj Seurat object
#' @param annotation_file Path to annotation RDS
#' @return Annotated and scored Seurat object
annotate_and_score <- function(chromatin_obj, annotation_file) {
  cat("Filtering to standard chromosomes...\n")
  grange_counts <- StringToGRanges(rownames(chromatin_obj), sep = c(":", "-"))
  grange_use <- seqnames(grange_counts) %in% standardChromosomes(grange_counts)
  chromatin_obj <- chromatin_obj[as.vector(grange_use), ]
  cat("Adding annotations...\n")
  annotations <- readRDS(annotation_file)
  Annotation(chromatin_obj) <- annotations
  genome(chromatin_obj) <- seqinfo(TxDb.Hsapiens.UCSC.hg38.knownGene)
  cat("Calculating TSS enrichment...\n")
  chromatin_obj <- TSSEnrichment(chromatin_obj)
  cat("Calculating nucleosome signal...\n")
  chromatin_obj <- NucleosomeSignal(chromatin_obj)
  return(chromatin_obj)
}

# ================================
# Main Script
# ================================

main <- function() {
  plan("multiprocess", workers = 4)
  options(future.globals.maxSize = 250000 * 1024^2)

  parser <- ArgumentParser(description = "Create Seurat Chromatin Object from binned genome intervals and fragment files.")
  parser$add_argument("--binnedgenomefilepath", type = "character", help = "Path to an rds/bed genomic interval file to get fragments per interval (bin)")
  parser$add_argument("--date", type = "character", help = "Date script is running; will be used in saving object names")
  parser$add_argument("--figdir", type = "character", help = "Figures directory; dir must already exist")
  parser$add_argument("--datadir", type = "character", help = "Data directory; dir must already exist")
  parser$add_argument("--workingdir", type = "character", help = "Working directory; dir must already exist")
  parser$add_argument("--makefragobj", type = "logical", help = "Whether to make fragment object (logical)")
  parser$add_argument("--fragmentfile", type = "character", help = "Fragment file path")
  parser$add_argument("--barcodesfile", type = "character", help = "Barcodes file path")
  parser$add_argument("--samplename", type = "character", help = "Sample name")
  parser$add_argument("--outname", type = "character", help = "Output name", required = FALSE)

  args <- parser$parse_args()
  cat(sprintf("Arguments: %s\n", paste(capture.output(str(args)), collapse = "; ")))

  setwd(args$workingdir)

  chromatin_object_base <- file.path(args$datadir, paste0("chromatin_objects-", args$date))
  annotation_file <- "~/common/annotations_UCSC_hg38_21092021.rds"

  binned_genome <- read_binned_genome(args$binnedgenomefilepath)
  barcodes <- read_barcodes(args$barcodesfile)

  if (isTRUE(args$makefragobj)) {
    chromatin_obj <- create_chromatin_object(args$fragmentfile, barcodes, binned_genome)
    saveRDS(chromatin_obj, paste0(chromatin_object_base, "_", args$samplename, "_", args$date, ".rds"))
  } else {
    chromatin_obj <- readRDS(paste0(chromatin_object_base, "_", args$samplename, "_", args$date, ".rds"))
  }

  chromatin_obj <- annotate_and_score(chromatin_obj, annotation_file)

  # Output file naming
  if (is.null(args$outname)) {
    save_outname <- paste0(chromatin_object_base, "_scores_", args$samplename, "_", args$date, ".rds")
    cat(sprintf("No outname supplied, using: %s\n", save_outname))
  } else {
    save_outname <- args$outname
    cat(sprintf("Using supplied outname: %s\n", save_outname))
  }
  saveRDS(chromatin_obj, save_outname)
  cat(sprintf("Finished script. Output saved to: %s\n", save_outname))
}

# ================================
# Run main
# ================================

if (sys.nframe() == 0) {
  tryCatch({
    main()
  }, error = function(e) {
    cat(sprintf("Script failed: %s\n", e$message))
    stop(e)
  })
}