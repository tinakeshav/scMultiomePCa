# ============================================================================
# Intra-patient TF variance analysis
# This script performs PCA on intra-patient TF variance and generates plots.
# PCA is run per sample (on Epithelial populations), then results are combined.
# Figures mostly used for F3 of manuscript
# ============================================================================


# =====================
# Load Libraries
# =====================
suppressPackageStartupMessages({
  library(argparse)
  library(Signac)
  library(Seurat)
  library(PCAtools)
  library(ggplot2)
  library(ggtext)
  library(dplyr)
  library(patchwork)
  library(ComplexHeatmap)
  library(circlize)
})


# =====================
# Argument
# =====================
parse_args <- function() {
  parser <- ArgumentParser(description = 'Intra-patient TF variance analysis')
  parser$add_argument('--multiome_files', nargs = '+', type = "character", required = TRUE,
                      help = "List of multiome RDS files (one per patient)")
  parser$add_argument('--patient_names', nargs = '+', type = "character", required = TRUE,
                      help = "Names of patients (order must match multiome_files)")
  parser$add_argument('--metadata_files', nargs = '+', type = "character", required = TRUE,
                      help = "List of metadata CSV files (one per patient, order must match multiome_files)")
  parser$add_argument('--cluster_ident', type = "character", default =  "wsnn_res.0.4", required = FALSE,
                      help = "Set identity to cluster")
  parser$add_argument('--output_prefix', type = "character", required = TRUE,
                      help = "Prefix for output files (e.g., figures/analysis)")
  parser$add_argument('--date', type = "character", required = TRUE,
                      help = "Date string for output files")
  args <- parser$parse_args()
  return(args)
}

# =====================
# Functions
# =====================

#' Calculate median chromVAR deviation per group
#' @param seurat_obj Seurat object
#' @param group_colname Column name in meta.data for grouping
#' @param chromvar_motif_list List of motif names
#' @return Dataframe of median chromVAR deviations per group
get_median_chromvar <- function(seurat_obj, group_colname, chromvar_motif_list) {
  median_dfs <- lapply(unique(seurat_obj@meta.data[[group_colname]]), function(grp) {
    grp_indices <- seurat_obj@meta.data[[group_colname]] == grp
    cells <- rownames(seurat_obj@meta.data[grp_indices, ])
    chromvar_data <- seurat_obj$chromvar@data[chromvar_motif_list, cells, drop = FALSE]
    median_chromvar <- as.data.frame(apply(chromvar_data, 1, median, na.rm = TRUE))
    colnames(median_chromvar) <- grp
    return(median_chromvar)
  })
  return(do.call(cbind, median_dfs))
}

#' Process each sample: read data, compute median chromVAR, save results
#' @param multiome_files List of RDS file paths
#' @param cluster_ident Cluster identity column
#' @param metadata_files List of metadata CSVs
#' @param out_date Date string
#' @param output_dir Output directory
#' @return List of output CSV filenames
process_samples <- function(multiome_files, cluster_ident, metadata_files, out_date, output_dir = "figures/") {
  median_chromvar_files <- character()
  for (i in seq_along(multiome_files)) {
    cat(sprintf('Processing %s...\n', multiome_files[i]))
    sample <- tryCatch(readRDS(multiome_files[i]), error = function(e) {
      cat(sprintf('ERROR: Failed to read RDS: %s\n', e$message)); stop(e)
    })
    mt <- tryCatch(read.csv(metadata_files[i], row.names = 1), error = function(e) {
      cat(sprintf('ERROR: Failed to read metadata: %s\n', e$message)); stop(e)
    })
    if (!all(rownames(mt) == colnames(sample))) {
      cat('ERROR: Metadata rownames do not match sample colnames.\n')
      stop('Metadata rownames do not match sample colnames.')
    }
    if (!all(rownames(mt) == rownames(sample@meta.data))) {
      cat('ERROR: Metadata rownames do not match sample meta.data rownames.\n')
      stop('Metadata rownames do not match sample meta.data rownames.')
    }
    sample@meta.data <- mt
    chromvar_tfs <- rownames(sample$chromvar@data)
    median_chromvar_dev <- get_median_chromvar(sample, cluster_ident, chromvar_tfs)
    filename <- sub('.rds$', paste0('_chromvar_Markers_', cluster_ident, '_Median_', out_date, '.csv'), multiome_files[i])
    write.csv(median_chromvar_dev, filename)
    median_chromvar_files <- c(median_chromvar_files, filename)
    rm(sample); gc()
  }
  return(median_chromvar_files)
}

#' Run PCA on Malignant Epithelial clusters and save results
#' @param median_chromvar_dev Median chromVAR matrix
#' @param motif_names Motif annotation dataframe
#' @param output_dir Output directory
#' @param sample_name Sample name
#' @param out_date Date string
#' @return PC loadings dataframe
run_malignant_pca <- function(median_chromvar_dev, motif_names, output_dir, sample_name, out_date) {
  malignant_cols <- grep('Malignant Epithelial', colnames(median_chromvar_dev), value = TRUE)
  if (length(malignant_cols) > 1) {
    median_malig <- median_chromvar_dev[, malignant_cols, drop = FALSE]
    median_malig_nona <- median_malig[, colSums(is.na(median_malig)) == 0, drop = FALSE]
    median_malig_nona_t <- t(median_malig_nona)
    annot_df <- data.frame(cluster = rownames(median_malig_nona_t))
    rownames(annot_df) <- rownames(median_malig_nona_t)
    p <- PCAtools::pca(median_malig_nona_t, metadata = annot_df)
    p$loadings$family <- motif_names[rownames(p$loadings), 'family']
    p$loadings$colour <- motif_names[rownames(p$loadings), 'colour']
    pc_loadings_filename <- file.path(output_dir, paste0(sample_name, '_ChromVar_PC_Loadings_', out_date, '.csv'))
    write.csv(p$loadings, pc_loadings_filename)
    var_explained <- round(p$variance, 1)
    pca_plot <- PCAtools::biplot(p, lab = rownames(annot_df), colby = NULL) +
      labs(x = paste0("PC1 (", var_explained[1], "%)"),
           y = paste0("PC2 (", var_explained[2], "%)")) +
      theme_classic()
    pca_plot_filename <- file.path(output_dir, paste0(sample_name, '_PCA_TF_', out_date, '.pdf'))
    pdf(pca_plot_filename, height = 5, width = 5)
    print(pca_plot)
    dev.off()
    return(p$loadings)
  } else {
    cat(sprintf('WARNING: No Malignant Epithelial clusters found for %s\n', sample_name))
    return(NULL)
  }
}

#' Normalize PC loadings (absolute, per PC)
get_normalized_loadings <- function(loadings) {
  abs_loadings <- abs(loadings)
  return(abs_loadings / sum(abs_loadings))
}

#' Get top features per PC by percentile
get_top_features <- function(normalized_loadings, percentile = 0.1) {
  thresholds <- apply(normalized_loadings, 2, function(x) quantile(x, 1 - percentile))
  top_features <- lapply(seq_len(ncol(normalized_loadings)), function(i) {
    col_data <- normalized_loadings[, i, drop = FALSE]
    rownames(col_data)[col_data >= thresholds[[i]]]
  })
  names(top_features) <- colnames(normalized_loadings)
  return(top_features)
}

#' Read PC loadings matrices for all patients
read_pc_matrices <- function(pc_loadings_csv, patient_names, top_pc = NA) {
  loadings_data <- list()
  for (i in seq_along(pc_loadings_csv)) {
    pc_loadings_data <- read.csv(pc_loadings_csv[i], row.names = 1)
    if (!is.na(top_pc)) {
      pc_loadings_data <- pc_loadings_data[, top_pc, drop = FALSE]
    }
    colnames(pc_loadings_data) <- paste0(colnames(pc_loadings_data), "_", patient_names[i])
    loadings_data[[i]] <- pc_loadings_data
  }
  names(loadings_data) <- patient_names
  return(loadings_data)
}

#' Combine PC matrices (intersection of TFs)
combine_pc_matrices <- function(pc_loadings_list) {
  rownames_to_keep <- Reduce(intersect, lapply(pc_loadings_list, rownames))
  combined_loadings <- do.call(cbind, lapply(pc_loadings_list, function(df) {
    df[rownames_to_keep, , drop = FALSE]
  }))
  return(combined_loadings)
}

#' Create heatmap of PC loadings
create_pc_loadings_heatmap <- function(heatmap_matrix, patient_names, motif_names_df) {
  motif_families <- motif_names_df[rownames(heatmap_matrix), "family", drop = FALSE]
  family_colors <- setNames(motif_names_df$colour, motif_names_df$family)
  row_annot <- HeatmapAnnotation(
    family = as.vector(motif_families[, 'family']),
    col = list(family = family_colors),
    show_legend = TRUE,
    which = "row"
  )
  color_palette <- colorRamp2(c(min(heatmap_matrix, na.rm = TRUE), 0, max(heatmap_matrix, na.rm = TRUE)),
                              c("blue", "white", "red"))
  ht <- Heatmap(
    matrix = heatmap_matrix,
    name = "PC Loading",
    show_row_names = TRUE,
    show_column_names = TRUE,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    left_annotation = row_annot,
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 8),
    width = unit(8, "inches"),
    height = unit(20, "inches"),
    col = color_palette
  )
  return(ht)
}

#' Get column annotation dataframe for heatmap
get_column_annot_df <- function(heatmap_matrix) {
  row_annot_cols_filename <- 'data/colours_23032023_with_samples_averaged_vals_03032023_CNVrep.csv'
  row_annot_col_df <- read.csv(row_annot_cols_filename)
  row_annot_list <- setNames(row_annot_col_df$colour, row_annot_col_df$annotation)
  pat_names_per_pc <- sapply(colnames(heatmap_matrix), function(x) unique(sub("\\..*", "", x)))
  col_annot <- data.frame(pat_names_per_pc, row.names = colnames(heatmap_matrix))
  col_annot$colour <- row_annot_list[col_annot$pat_names_per_pc]
  return(col_annot)
}

#' Create row-sorted heatmap of PC loadings
create_pc_loadings_row_sorted_heatmap <- function(heatmap_matrix, patient_names, motif_names_df, tf_row_order) {
  motif_families <- motif_names_df[tf_row_order, "family", drop = FALSE]
  family_colors <- setNames(motif_names_df$colour, motif_names_df$family)
  row_annot <- HeatmapAnnotation(
    family = as.vector(motif_families[, 'family']),
    col = list(family = family_colors),
    show_legend = TRUE,
    which = "row"
  )
  col_annot_df <- get_column_annot_df(heatmap_matrix)
  col_annot <- HeatmapAnnotation(
    patient = as.vector(col_annot_df[, 'pat_names_per_pc']),
    col = list(patient = setNames(col_annot_df$colour, col_annot_df$pat_names_per_pc)),
    show_legend = TRUE,
    which = "col"
  )
  color_palette <- colorRamp2(c(0, 0.001, 0.01), c("#030389", "white", "#800404"))
  ht <- Heatmap(
    matrix = heatmap_matrix[tf_row_order, ],
    column_split = col_annot_df$pat_names_per_pc,
    name = "TF weight",
    show_row_names = TRUE,
    show_column_names = TRUE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    left_annotation = row_annot,
    top_annotation = col_annot,
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 14),
    column_names_gp = gpar(fontsize = 14),
    col = color_palette
  )
  return(ht)
}

#' Calculate total variance explained by families in top hits
get_family_variance <- function(normalized_loadings, motif_names_df) {
  tf_families <- motif_names_df[rownames(normalized_loadings), "family"]
  family_variance <- aggregate(normalized_loadings, by = list(Family = tf_families), FUN = sum, na.rm = TRUE)
  family_matrix <- as.matrix(family_variance[, -1])
  rownames(family_matrix) <- family_variance$Family
  return(family_matrix)
}

#' Create heatmap for family variance
create_family_variance_heatmap <- function(family_matrix, motif_names) {
  unique_families <- unique(motif_names$family)
  family_colors <- setNames(unique(motif_names$colour[match(unique_families, motif_names$family)]), unique_families)
  row_annot <- HeatmapAnnotation(
    family = rownames(family_matrix),
    col = list(family = family_colors),
    show_legend = FALSE,
    which = "row"
  )
  col_annot_df <- get_column_annot_df(family_matrix)
  col_annot <- HeatmapAnnotation(
    patient = as.vector(col_annot_df[, 'pat_names_per_pc']),
    col = list(patient = setNames(col_annot_df$colour, col_annot_df$pat_names_per_pc)),
    show_legend = FALSE,
    which = "col"
  )
  color_palette <- colorRamp2(c(0, 0.001, max(family_matrix)), c("#030389", "white", "#800404"))
  ht <- Heatmap(
    matrix = family_matrix,
    column_split = col_annot_df$pat_names_per_pc,
    name = "Variance Explained",
    show_row_names = TRUE,
    show_column_names = TRUE,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    left_annotation = row_annot,
    top_annotation = col_annot,
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 14),
    column_names_gp = gpar(fontsize = 14),
    col = color_palette
  )
  return(ht)
}

#' Plot TF family composition in PC bins
plot_pc_bin_tf_families <- function(loadings_df, pc_column = "PC1", bin_number = 1, n_bins = 30) {
  pc_loadings <- loadings_df %>%
    arrange(desc(!!sym(pc_column))) %>%
    select(c(!!sym(pc_column), family, colour)) %>%
    mutate(bin = ceiling(row_number() / (n() / n_bins)))
  bin_data <- pc_loadings %>%
    filter(bin == bin_number) %>%
    arrange(!!sym(pc_column)) %>%
    mutate(num = seq(1, n())) %>%
    mutate(num = factor(num, levels = unique(num)))
  bin_data$tf_name <- rownames(bin_data)
  tf_family_colours <- setNames(bin_data$colour, bin_data$family)
  bar_plot <- ggplot(bin_data) +
    geom_col(aes(x = 1, y = num, fill = family), position = "stack", width = 1) +
    geom_text(aes(x = 1, y = num, label = tf_name), hjust = 0, size = 3) +
    scale_fill_manual(values = tf_family_colours) +
    theme_void() +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
          axis.ticks = element_blank(), axis.text = element_blank()) +
    coord_cartesian(clip = 'off', xlim = c(0, 3)) +
    labs(title = sprintf("%s Bin %d/%d", pc_column, bin_number, n_bins))
  return(bar_plot)
}

# =====================
# Run main
# =====================
main <- function() {
  args <- parse_args()
  multiome_files <- args$multiome_files
  patient_names <- args$patient_names
  metadata_files <- args$metadata_files
  cluster_ident <- args$cluster_ident
  out_date <- args$date
  output_prefix <- args$output_prefix

  motif_name_filename <- 'data/JASPAR2020_tag_motif_info_colours_edited.csv'
  motif_names <- read.csv(motif_name_filename, row.names = 1)

  median_chromvar_files <- process_samples(multiome_files, cluster_ident, metadata_files, out_date, output_dir = "figures/")
  pc_loadings_list <- list()
  for (i in seq_along(median_chromvar_files)) {
    sample_name <- patient_names[i]
    pc_loadings <- run_malignant_pca(median_chromvar_dev, motif_names, output_dir = "figures/", sample_name = sample_name, out_date = out_date)
    if (!is.null(pc_loadings)) {
      pc_loadings_list[[sample_name]] <- pc_loadings
    }
  }
  names(pc_loadings_list) <- patient_names

  # Per-sample PC bin TF family plots
  h <- 10; w <- 10
  for (i in seq_along(pc_loadings_list)) {
    loadings_df <- pc_loadings_list[[i]]
    plot_first_pc1 <- plot_pc_bin_tf_families(loadings_df, "PC1", 1)
    plot_last_pc1 <- plot_pc_bin_tf_families(loadings_df, "PC1", 30)
    filename_first <- file.path("figures", paste0(sample_name, '_PCA_TF_gradients_PC1_first_bin_', out_date, '.pdf'))
    filename_last <- file.path("figures", paste0(sample_name, '_PCA_TF_gradients_PC1_last_bin_', out_date, '.pdf'))
    pdf(filename_first, height = h, width = w); print(plot_first_pc1); dev.off()
    pdf(filename_last, height = h, width = w); print(plot_last_pc1); dev.off()
    plot_first_pc2 <- plot_pc_bin_tf_families(loadings_df, "PC2", 1)
    plot_last_pc2 <- plot_pc_bin_tf_families(loadings_df, "PC2", 30)
    filename_first_pc2 <- file.path("figures", paste0(sample_name, '_PCA_TF_gradients_PC2_first_bin_', out_date, '.pdf'))
    filename_last_pc2 <- file.path("figures", paste0(sample_name, '_PCA_TF_gradients_PC2_last_bin_', out_date, '.pdf'))
    pdf(filename_first_pc2, height = h, width = w); print(plot_first_pc2); dev.off()
    pdf(filename_last_pc2, height = h, width = w); print(plot_last_pc2); dev.off()
    combined_plot <- (plot_first_pc1 + plot_last_pc1) / (plot_first_pc2 + plot_last_pc2)
    filename_combined <- file.path("figures", paste0(sample_name, '_PCA_TF_gradients_combined_', out_date, '.pdf'))
    pdf(filename_combined, height = h * 2, width = w * 2); print(combined_plot); dev.off()
  }

  #combine loadings for gloal fig
  combined_loadings <- as.data.frame(combine_pc_matrices(pc_loadings_list))
  normalized_loadings <- as.data.frame(apply(combined_loadings, 2, get_normalized_loadings))

  #global heatmap
  pc_heatmap <- create_pc_loadings_heatmap(normalized_loadings, patient_names, motif_names)
  pdf(paste0("figures/pc_loadings_heatmap_", out_date, ".pdf"), width = 20, height = 30)
  print(pc_heatmap)
  dev.off()

  # Top features
  top_tf_per_pc <- get_top_features(normalized_loadings, 0.2)
  top_tf_per_pc_df <- data.frame(PC = names(unlist(top_tf_per_pc)), Features = unlist(top_tf_per_pc))
  top_tf_per_pc_df$pat <- sub("\\..*", "", top_tf_per_pc_df$PC)
  top_tf_per_pc_df_agg <- aggregate(PC ~ Features, data = top_tf_per_pc_df, FUN = function(x) paste(x, collapse = ","))
  top_tf_per_pc_df_agg$feature_repeat <- sapply(top_tf_per_pc_df_agg$PC, function(x) length(strsplit(x, ',')[[1]]))
  splitted <- sapply(top_tf_per_pc_df_agg$PC, function(x) strsplit(x, ',')[[1]])
  top_tf_per_pc_df_agg$patients <- sapply(splitted, function(x) unique(sub("\\..*", "", x)))
  top_tf_per_pc_df_agg$patient_repeat <- sapply(top_tf_per_pc_df_agg$patients, function(x) length(x))
  top_tf <- unique(unlist(top_tf_per_pc))

  pc_heatmap <- create_pc_loadings_heatmap(normalized_loadings[top_tf, ], patient_names, motif_names)
  pdf(paste0("figures/pc_top_loadings_heatmap_", out_date, ".pdf"), width = 20, height = 30)
  print(pc_heatmap)
  dev.off()

  # sort top features by their family. 
  # i compiled family from the JASPAR2020 database
  top_tf_per_pc_df_agg <- top_tf_per_pc_df_agg %>% mutate(family = motif_names[Features, 'family']) %>% arrange(desc(feature_repeat), family)
  family_order <- unique(top_tf_per_pc_df_agg$family)
  top_tf_per_pc_df_agg$family <- factor(top_tf_per_pc_df_agg$family, levels = family_order)
  top_tf_per_pc_df_agg <- top_tf_per_pc_df_agg %>% arrange(family, desc(feature_repeat))

  pc_heatmap <- create_pc_loadings_row_sorted_heatmap(normalized_loadings[top_tf_per_pc_df_agg$Features, ], patient_names, motif_names, tf_row_order = top_tf_per_pc_df_agg$Features)
  pdf(paste0("figures/pc_top_loadings_heatmap_sorted_", out_date, ".pdf"), width = 20, height = 30)
  print(pc_heatmap)
  dev.off()

  # Family variance
  family_variance_df <- get_family_variance(normalized_loadings[top_tf, ], motif_names)
  family_variance_heatmap <- create_family_variance_heatmap(family_variance_df, motif_names)
  pdf(paste0("figures/pc_family_variance_heatmap_", out_date, ".pdf"), width = 12, height = 16)
  print(family_variance_heatmap)
  dev.off()
  write.csv(family_variance_df, file = paste0("data/pc_family_variance_", out_date, ".csv"))

  #  Plot selected TFs of interest
  tfs_of_interest <- c('FOS', 'FOSL1', 'BACH1', 'BATF', 'STAT1::STAT2', 'FOXP2', 'FOXA2', 'FOXA1', 'GATA3', 'GATA5',
                      'HOXC13', 'HOXA13', 'IRF1', 'NR3C1', 'MAX', 'MYC', 'SNAI1', 'NFYA', 'NKX6-1', 'NFIC', 'ARNT2',
                      'TEAD3', 'TCF4', 'GRHL1', 'KLF5', 'CEBPA', 'ASCL1(var.2)', 'PAX6', 'RUNX2', 'TWIST1', 'SP4',
                      'YY1', 'HIF1A')
  top_tf_per_pc_df_agg_select <- top_tf_per_pc_df_agg[tfs_of_interest, ] %>% arrange(family, desc(feature_repeat))
  pc_heatmap <- create_pc_loadings_row_sorted_heatmap(normalized_loadings[top_tf_per_pc_df_agg_select$Features, ],
                                                      patient_names, motif_names, tf_row_order = top_tf_per_pc_df_agg_select$Features)
  pdf(paste0("figures/pc_top_loadings_heatmap_sorted_select_tfs_", out_date, ".pdf"), width = 10, height = 15)
  print(pc_heatmap)
  dev.off()
}

# =====================
# Run main if script is executed with arguments
# =====================
if (sys.nframe() == 0) {
  main()
}