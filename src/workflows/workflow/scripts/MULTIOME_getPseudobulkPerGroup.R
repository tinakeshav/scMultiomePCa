# ================================
# Description:
# Pseudobulk sample counts for groups for peaks called by intrapat
# ================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(DESeq2)
  library(argparse)
})


#' Print log messages with timestamp
log_message <- function(...) {
  cat(sprintf("[%s] ", Sys.time()), ..., "\n")
}

#' Get pseudobulked raw counts per group for a Seurat object
#'
#' This function aggregates raw counts for each group defined by a metadata column in a Seurat object.
#'
#' @param s Seurat object
#' @param group_colname Character. Metadata column name to group by.
#' @param assay Character. Assay name to extract counts from.
#' @return Data frame of pseudobulked raw counts (features x groups)
get_pseudobulk_counts = function(s, group_colname, assay) {
	# returns a dataframe of pseudobulked raw counts per assay
	# <s> : Seurat object
	# <group_colname> : metadata column name to subset by
	pseudo_dfs = list()

	# gets unique values per group
	for(grp in unique(s@meta.data[, group_colname])){
	    grp_indices = s@meta.data[, group_colname] == grp 
	    cells = rownames(s@meta.data[grp_indices, ])

	    group_counts = s[[assay]]@counts[, cells]
	    group_pseudobulk = as.data.frame(Matrix::rowSums(group_counts))
	    colnames(group_pseudobulk) = grp
	    pseudo_dfs = c(pseudo_dfs, list(group_pseudobulk))
	}
	return(do.call(cbind, pseudo_dfs))
}

#' Normalize pseudobulked counts using DESeq2 variance stabilizing transformation
#'
#' This function takes a pseudobulked count data frame and returns a DESeq2 object
#' with variance-stabilized counts.
#'
#' @param pseudobulk_counts_df Data frame of pseudobulked counts (features x groups)
#' @return DESeq2 object with VST-normalized counts
#' @importFrom DESeq2 DESeqDataSetFromMatrix vst
#' @export
deseq_norm_pseudobulk_counts = function(pseudobulk_counts_df){
	coldata = data.frame(names = names(pseudobulk_counts_df), 
		     row.names = names(pseudobulk_counts_df))
	dds <- DESeqDataSetFromMatrix(countData = pseudobulk_counts_df,
		      colData = coldata,
		      design = ~1)
	vst_dds = vst(dds)
	return(vst_dds)
}

#' Subset normalized counts to a set of peaks from a BED file
#'
#' This function reads a BED file and subsets the normalized count matrix to only those peaks.
#'
#' @param peak_bedfilename Path to BED file with peaks
#' @param vst_counts Matrix of normalized counts (features x groups)
#' @return Subsetted matrix of normalized counts (peaks x groups)
subsample_on_peaks = function(peak_bedfilename, vst_counts){
	peaks = read.table(peak_bedfilename)
	peak_featurenames = paste(peaks[,'V1'], peaks[,'V2']+1, peaks[,'V3'], sep = '-')
	return(vst_counts[peak_featurenames, ])
}

#' Load Seurat object and (optionally) metadata
#'
#' Loads a Seurat object from an RDS file and, if provided, replaces its metadata with a CSV file.
#' Checks for file existence and matching row/column names.
#'
#' @param sample_filename Path to Seurat RDS file
#' @param sample_metadata_filename Path to metadata CSV file (optional)
#' @return Seurat object with updated metadata (if provided)
load_sample <- function(sample_filename, sample_metadata_filename = NULL) {
  if (!file.exists(sample_filename)) {
    stop(paste('Sample file does not exist:', sample_filename))
  }
  sample <- readRDS(sample_filename)
  if (!is.null(sample_metadata_filename) && sample_metadata_filename != '') {
    if (!file.exists(sample_metadata_filename)) {
      stop(paste('Sample metadata file does not exist:', sample_metadata_filename))
    }
    sample_metadata <- read.csv(sample_metadata_filename, row.names = 1, header = TRUE, check.names = FALSE)
    if (!all(rownames(sample_metadata) == colnames(sample))) {
      stop('Row names of metadata do not match column names of Seurat object.')
    }
    sample@meta.data <- sample_metadata
  }
  return(sample)
}

#' Subset Seurat object by metadata column and value
#'
#' Subsets a Seurat object to only those cells where a metadata column matches a given value.
#'
#' @param sample Seurat object
#' @param col_to_filteron Character. Metadata column to filter on.
#' @param col_to_filteron_value Character. Value to filter for.
#' @return Subsetted Seurat object
subset_sample <- function(sample, col_to_filteron, col_to_filteron_value) {
  if (!is.null(col_to_filteron) && col_to_filteron != '' &&
      !is.null(col_to_filteron_value) && col_to_filteron_value != '') {
    if (!(col_to_filteron %in% colnames(sample@meta.data))) {
      stop(paste('Column to filter on not found in metadata:', col_to_filteron))
    }
    log_message('Subsetting sample on', col_to_filteron, '=', col_to_filteron_value)
    cells <- sample@meta.data[, col_to_filteron] == col_to_filteron_value
    sample <- subset(sample, cells = rownames(sample@meta.data[cells, ]))
    log_message('Sample dimensions after subsetting:', paste(dim(sample), collapse = ' x '))
  }
  return(sample)
}

#' Process a single assay: pseudobulk, normalize, correlate, and save results
#'
#' For a given assay, this function computes pseudobulked counts, normalizes them,
#' computes correlation matrices, and writes all outputs to disk. Optionally, it
#' recomputes correlations for promoter/non-promoter peaks if provided.
#'
#' @param sample Seurat object
#' @param colname Character. Metadata column to group by.
#' @param assay Character. Assay name.
#' @param pseudobulk_countmat_prefix Character. Prefix for output files.
#' @param date Character. Date string for filenames.
#' @param promoter_peaks Character. Path to promoter peaks BED file (optional).
#' @param non_promoter_peaks Character. Path to non-promoter peaks BED file (optional).
#' @return None. Writes files to disk.
process_assay <- function(sample, colname, assay, pseudobulk_countmat_prefix, date, promoter_peaks, non_promoter_peaks) {
  pseudo_df <- get_pseudobulk_counts(sample, colname, assay)
  fname_unnorm <- paste0(pseudobulk_countmat_prefix, '_', assay, '_unnormPseudoCounts_', date, '.csv')
  log_message('Saving unnormalized pseudobulk counts to', fname_unnorm)
  write.csv(pseudo_df, fname_unnorm)

  vst_dds <- deseq_norm_pseudobulk_counts(pseudo_df)
  fname_vst_rds <- paste0(pseudobulk_countmat_prefix, '_', assay, '_DESeq2vstPseudoCounts_', date, '.rds')
  fname_vst_csv <- paste0(pseudobulk_countmat_prefix, '_', assay, '_DESeq2vstPseudoCounts_', date, '.csv')
  log_message('Saving DESeq2 VST normalized counts to', fname_vst_rds, 'and', fname_vst_csv)
  saveRDS(vst_dds, fname_vst_rds)
  write.csv(assay(vst_dds), fname_vst_csv)

  cor_mat <- cor(assay(vst_dds))
  fname_cor <- paste0(pseudobulk_countmat_prefix, '_', assay, '_DESeq2vstCorrMat_', date, '.csv')
  log_message('Saving correlation matrix to', fname_cor)
  write.csv(cor_mat, fname_cor)

  if (!is.null(promoter_peaks) && promoter_peaks != '') {
    subset_peaks_norm <- subsample_on_peaks(promoter_peaks, assay(vst_dds))
    cor_mat <- cor(subset_peaks_norm)
    fname_cor_prom <- paste0(pseudobulk_countmat_prefix, '_', assay, '_DESeq2vstCorrMat_PromoterPeaks', date, '.csv')
    log_message('Saving promoter peaks correlation matrix to', fname_cor_prom)
    write.csv(cor_mat, fname_cor_prom)
  }
  if (!is.null(non_promoter_peaks) && non_promoter_peaks != '') {
    subset_peaks_norm <- subsample_on_peaks(non_promoter_peaks, assay(vst_dds))
    cor_mat <- cor(subset_peaks_norm)
    fname_cor_nonprom <- paste0(pseudobulk_countmat_prefix, '_', assay, '_DESeq2vstCorrMat_NonPromoterPeaks', date, '.csv')
    log_message('Saving non-promoter peaks correlation matrix to', fname_cor_nonprom)
    write.csv(cor_mat, fname_cor_nonprom)
  }
}


main <- function() {
  options(scipen = 999)
  parser <- ArgumentParser(description = 'Pseudobulk sample counts for groups per assay')
  parser$add_argument('--sample_filename', type = 'character', required = TRUE, help = 'Seurat object RDS file')
  parser$add_argument('--sample_metadata_filename', type = 'character', required = FALSE, help = 'Sample metadata CSV file')
  parser$add_argument('--colname', type = 'character', required = TRUE, help = 'Column name to group by and save pseudobulked counts')
  parser$add_argument('--assays', type = 'character', nargs = '+', required = TRUE, help = 'Assays to get counts for and normalize')
  parser$add_argument('--pseudobulk_countmat_prefix', type = 'character', required = TRUE, help = 'Prefix for pseudobulked count matrix files')
  parser$add_argument('--date', type = 'character', required = TRUE, help = 'Date to add before suffix when saving')
  parser$add_argument('--col_to_filteron', type = 'character', required = FALSE, help = 'Column to filter on prior to pseudobulk')
  parser$add_argument('--col_to_filteron_value', type = 'character', required = FALSE, help = 'Only retrieve these values')
  parser$add_argument('--promoter_peaks', type = 'character', required = FALSE, help = 'Promoter peak regions to filter')
  parser$add_argument('--non_promoter_peaks', type = 'character', required = FALSE, help = 'Non-promoter peak regions to filter')
  parser$add_argument('--output_txtfile', type = 'character', required = FALSE, help = 'Output txt file for Snakemake pipeline')
  args <- parser$parse_args()

  log_message('Parsed arguments:')
  print(args)

  if (length(args$assays) > 1 && (!is.null(args$promoter_peaks) && args$promoter_peaks != '')) {
    stop('Can only perform subsetting on a single assay at a time when promoter_peaks is supplied.')
  }

  sample <- load_sample(args$sample_filename, args$sample_metadata_filename)
  sample <- subset_sample(sample, args$col_to_filteron, args$col_to_filteron_value)

  for (assay in args$assays) {
    process_assay(
      sample = sample,
      colname = args$colname,
      assay = assay,
      pseudobulk_countmat_prefix = args$pseudobulk_countmat_prefix,
      date = args$date,
      promoter_peaks = args$promoter_peaks,
      non_promoter_peaks = args$non_promoter_peaks
    )
  }

  if (!is.null(args$output_txtfile) && args$output_txtfile != '') {
    log_message('Writing success stamp to', args$output_txtfile)
    writeLines(args$pseudobulk_countmat_prefix, args$output_txtfile)
  }
  log_message('Pseudobulk processing completed successfully.')
}

# -----------------------------
# Run main if script is executed with command line arguments
# -----------------------------

if (sys.nframe() == 0) {
  tryCatch(
    main(),
    error = function(e) {
      cat(sprintf("[ERROR %s] %s\n", Sys.time(), e$message))
      quit(status = 1)
    }
  )
}
# Results can be used to plot promoter and non-promoter heatmaps!