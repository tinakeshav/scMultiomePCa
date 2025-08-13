# Heatmaps for populations
# Functional script currently designed for figure 1
# Can be modified to be used for figure 3

#' Get correlation filenames
#' @param pseudobulk_countmat_prefix prefix for pseudobulk count matrices, used to obtain output correlation filenames
#' @param date date for filenames that get saved, must be the same date as the pseudobulk count matrices
#' @return list of correlation filenames
get_corr_filenames <- function(pseudobulk_countmat_prefix, date) {
    fname_cor_prom <- paste0(pseudobulk_countmat_prefix, '_', 'intrapat_peaks', '_DESeq2vstCorrMat_PromoterPeaks', date, '.csv')
    fname_cor_nonprom <- paste0(pseudobulk_countmat_prefix, '_', 'intrapat_peaks', '_DESeq2vstCorrMat_NonPromoterPeaks', date, '.csv')
    fname_vst_csv <- paste0(pseudobulk_countmat_prefix, '_', 'RNA', '_DESeq2vstPseudoCounts_', date, '.csv')
	return(assay_corr_filenames)
}

#' Get counts for each population
#' @param sample_metadata_filename filename for sample metadata
#' @param colname column name for the population to get counts for; MUST be the same as column name in correlation files or it will break.
#' @return counts for each population
get_counts <- function(sample_metadata_filename, colname) {
	mt = read.csv(sample_metadata_filename)
	counts = as.data.frame(table(mt$colname)) %>% rename(count = Freq)
	rownames(counts) = counts$Var1
	return(counts)
}
###############################################################################
# Load Arguments and Libraries
###############################################################################
library(argparse)

parser <- ArgumentParser(description='Plotting phenotypic corr') 


parser$add_argument('--pseudobulk_countmat_prefix', type = "character",
				required = TRUE, help = "prefix for pseudobulk count matrices, used to obtain output correlation filenames")

parser$add_argument('--date', type = "character",
		    required = TRUE, help = "date for filenames that get saved, must be the same date as the pseudobulk count matrices")

parser$add_argument('--prefix_filename_for_figures', type = "character",
		    required = TRUE, help = "prefix for saving figures")

parser$add_argument('--sample_metadata_filename', type = "character",
				required = TRUE, help = "filename for sample metadata, used to obtain count matrix for right most column of plot")

parser$add_argument('--row_annot_cols_filename', type = "character", 
		    default = 'data/colours_23032023_with_samples_averaged_vals_03032023_CNVrep.csv',
		    required = FALSE, help = "filename for annotations")
parser$add_argument('--colname', type = "character",
		    required = TRUE, help = "column name for the population to get counts for; MUST be the same as column name in correlation files or it will break.")
parser$add_argument('--output_txtfile', type = "character",
				required = TRUE, help = "filename for output text file that indicates completion of script to snakemake")

args <- parser$parse_args()


date = args$date
prefix_filename_for_figures = args$prefix_filename_for_figures 
sample_metadata_filename = args$sample_metadata_filename
row_annot_cols_filename = args$row_annot_cols_filename
colname = args$colname


# Load libraries
library(ComplexHeatmap)
library(circlize)

###############################################################################
# Read matrices
###############################################################################
# Read first matrix to plot. 

assay_corr_filenames = get_corr_filenames(args$pseudobulk_countmat_prefix)

all_dfs = lapply(assay_corr_filenames, 
	         read.csv, 
	         row.names = 1, 
	         header = TRUE, 
	         check.names = FALSE)

corr_mat = read.csv(assay_corr_filenames[[1]],
	            row.names = 1, 
	            header = TRUE, 
	            check.names = FALSE)
print(rownames(corr_mat))

# Plot first heatmap and get row order
ht_init = Heatmap(corr_mat)
rowname_order = rownames(corr_mat)[row_order(ht_init)]
print(rowname_order)

# Rename columns and rows in our matrix
# Get row annot colours
row_annot_col_df = read.csv(row_annot_cols_filename)
row_annot_list = setNames(row_annot_col_df$colour, 
			  row_annot_col_df$annotation)

# Get population num and annot
population_num = sub("_.*$", "", rowname_order)
population_annot = sub("^[^_]*_", "", rowname_order)
population_annot_cols = row_annot_list[population_annot]
print(population_annot_cols)
# Retrieve 
annot_df = data.frame(population_num = population_num, 
		       population_annot = population_annot,
		       row.names = rowname_order)

counts_info = get_counts(sample_metadata_filename, colname)
annot_df$count = counts_info[rownames(annot_df), 'count']
print(annot_df)


###############################################################################
# Plot the first set of heatmaps
###############################################################################
size_div = 1.8

# non prom atac colours
assay_diag_col = colorRamp2(c(min(all_dfs[[1]]), max(all_dfs[[1]])),
			    c('#0464d0', '#ffd434'))
print(min(all_dfs[[1]]))
print(max(all_dfs[[1]]))
print(min(all_dfs[[2]]))
print(max(all_dfs[[2]]))
print(min(all_dfs[[3]]))
print(max(all_dfs[[3]]))
# prom atac colours
prom_ATAC_col = colorRamp2(c(min(all_dfs[[2]]), max(all_dfs[[2]])),
		           c('#0464d0', '#ffd434'))
# taken from TCGA
RNA_col = colorRamp2(c(min(all_dfs[[3]]), max(all_dfs[[3]])),
		    c('#293398', '#f8ef13'))
 

row_a = HeatmapAnnotation(
		   clusters = anno_text(annot_df[rowname_order,'population_num']) 
		   , ht_annot = as.vector(annot_df[rowname_order,'population_annot'])
		   , col = list(ht_annot = row_annot_list[annot_df[,'population_annot']])
		   , show_legend = c(ht_annot = FALSE)
		   , annotation_name_gp = gpar(fontsize = 0)
		   , which = 'row'
		   )
print('done ra')

col_a = HeatmapAnnotation(
		   ht_annot = as.vector(annot_df[rowname_order,'population_annot'])
		   , clusters = anno_text(annot_df[rowname_order,'population_num']) 
		   , col = list(ht_annot = row_annot_list[annot_df[,'population_annot']])
		   , show_legend = c(ht_annot = FALSE)
		   , annotation_name_gp = gpar(fontsize = 0)
		   , which = 'column'
		   )
print('done col')

row_a_last = HeatmapAnnotation(
		  Count = anno_barplot(annot_df[, 'count']
			 , bar_width = 0.9 
			 , border = TRUE
			 , width = unit(3, 'cm')
			 , axis_param = list(gp = gpar(fontsize = 12))
			 )
		   , which = 'row'
		   )


corr_ht = Heatmap(corr_mat[rowname_order, rowname_order]
		      , cluster_columns = FALSE
		      , cluster_rows = TRUE 
		      , col = assay_diag_col
		      , show_row_name = FALSE
		      , show_column_name = FALSE
		      #, show_heatmap_legend = FALSE 
		      , heatmap_legend_param = list(legend_direction = 'horizontal'
						    , title = 'Non-promoter ATAC'
						    , legend_height = unit(2, "cm")
						    , legend_width = unit(3, "cm")
						    )
		      , left_annotation = row_a
		      , bottom_annotation = col_a
		      , width = 9/size_div
		      , height = 10/size_div
		      ) 

prom_atac_ht = Heatmap(all_dfs[[2]][rowname_order, rowname_order]
		      , cluster_columns = FALSE
		      , cluster_rows = FALSE
		      , col = prom_ATAC_col
		      , show_row_name = FALSE
		      , show_column_name = FALSE
		      , bottom_annotation = col_a
		      #, show_heatmap_legend = FALSE 
		      , heatmap_legend_param = list(legend_direction = 'horizontal'
						    , title = 'Promoter ATAC'
						    , legend_height = unit(2, "cm")
						    , legend_width = unit(3, "cm")
						    )
		      , width = 9/size_div
		      , height = 10/size_div
		      ) 

rna_ht = Heatmap(all_dfs[[3]][rowname_order, rowname_order]
		      , cluster_columns = FALSE
		      , cluster_rows = FALSE
		      , col = RNA_col
		      , show_row_name = FALSE
		      , show_column_name = FALSE
		      , bottom_annotation = col_a
		      #, show_heatmap_legend = FALSE 
		      , heatmap_legend_param = list(legend_direction = 'horizontal'
						    , title = 'RNA'
						    , legend_height = unit(2, "cm")
						    , legend_width = unit(3, "cm")
						    )
		      , right_annotation = row_a_last
		      , width = 9/size_div
		      , height = 10/size_div
		      ) 

filename = paste0(prefix_filename_for_figures, '_nonPromPromRNA_CorrsHeatmaps_', date, '.pdf')
pdf(filename, height = 10/size_div, width = 35/size_div)
corr_ht + prom_atac_ht + rna_ht 
dev.off()

# output a text file if script runs without errors
if (!is.null(args$output_txtfile) && args$output_txtfile != '') {
  writeLines("success", args$output_txtfile)
}
