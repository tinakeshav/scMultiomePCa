# R script to plot data from scripts/clonal_js_distance.py
source('scripts/utils_plot_heatmaps.R')
library(ggplot2)

##############################################################
# Load Arguments
##############################################################
library(argparse)

parser <- ArgumentParser(description='heatmaps for phys/hamm distance; on name')

parser$add_argument('--date', type = "character",
		    required = TRUE, help = "date for filenames that get saved")

parser$add_argument('--prefix_filename_for_figures', type = "character",
		    required = TRUE, help = "prefix for saving figures")

parser$add_argument('--colours_df_filename', type = "character",
		    required = TRUE, help = "prefix for saving figures")

parser$add_argument('--js_distance_filename', type = "character", 
		    required = TRUE, help = "JS distance csv file")

parser$add_argument('--heatmaptype', type = "character", 
		    required = FALSE, help = "heatmap type for annotations, \\
		  			      see LBN cluster_annot_short, new_name_mapping, count")
		  			      
parser$add_argument('--csv_file_with_mapping_and_count', type = "character", 
		    required = FALSE, help = "csv file for metadata column w/ 3 columns, \\
		  			      see LBN cluster_annot_short, new_name_mapping, count")

parser$add_argument('--filter_name_on_pattern', type = "character"
		    , required = FALSE, help = "filters on name pattern")

parser$add_argument('--filter_name_on_pattern_absence'
		    , default = 'FALSE', type = "character"
		    , required = FALSE, help = "filters on opposite of \\
		    			     <filer_name_on_pattern> \\
					     can pass anything")

parser$add_argument('--filter_on_count', type = "double", 
		    required = FALSE, help = "filters on count")


args <- parser$parse_args()


date = args$date
prefix_filename_for_figures = args$prefix_filename_for_figures
colours_df_filename = args$colours_df_filename 
js_distance_filename = args$js_distance_filename
heatmaptype = args$heatmaptype
cluster_csv_file = args$csv_file_with_mapping_and_count
pattern_to_filter = args$filter_name_on_pattern
pattern_not_to_filter = args$filter_name_on_pattern_absence
count_to_filter = args$filter_on_count


################################################################
# Function for cleaner code 
################################################################
row_fontsize = 20 
col_fontsize = 20 
fsize = 18
################################################################
# Plotting JS distanc  
################################################################
# Lastly, plot clonal distributions as annotations
# Read clone by sample proportions
js_df = read.csv(js_distance_filename,
			   check.names = F,
			   row.names = 1,
			   header = TRUE)
print(js_df)


# Obtain colours for filling in annotations
colours_df = read.csv(colours_df_filename, check.names = F, header = TRUE)
colours_mapping = setNames(colours_df[,'colour'], colours_df[,'annotation'])

###############################################################################
# Plot the other matrix, if available
###############################################################################

# get clusters to plot
clusters_to_plot = cluster_csv_file
if(length(cluster_csv_file) > 0){
	clust_info = read.csv(cluster_csv_file, 
				row.names = 1,
				header = TRUE)
	if(length(pattern_to_filter) > 0){
		if(pattern_not_to_filter  == 'TRUE'){
			grepped_rows = grep(pattern_to_filter
					    , row.names(clust_info)
					    , invert = TRUE)
			clust_info = clust_info[grepped_rows,, drop = FALSE]
			print('subsetted cluster info on pattern')
			}
		else{
			grepped_rows = grep(pattern_to_filter, row.names(clust_info))
			clust_info = clust_info[grepped_rows,, drop = FALSE]
			print('subsetted cluster info on pattern')
		}
	}
	if(length(count_to_filter) > 0){
		clust_info = clust_info[clust_info[, 'cluster_annot_short'] > count_to_filter,, drop = FALSE]
		print('subsetted cluster info on count')
	}
	# reset value of clusters_to_plot
	clusters_to_plot = rownames(clust_info)
    print('clusters to plot')
    print(clusters_to_plot)
    js_df = js_df[clusters_to_plot, clusters_to_plot]
}


# Proportion colour mappings from fig_plot_population_heatmaps.R
col_mapping = colorRamp2(c(0,1), 
			c('#c5dae8', "#204f6f"))

ht_population_js = Heatmap(js_df
		       , show_row_names = FALSE 
		       , row_names_side = 'right'
		       , show_column_names = FALSE 
		       , cluster_columns = TRUE 
		       , col = col_mapping
		       , show_heatmap_legend = TRUE 
		       , width = 4
		       , row_names_gp = gpar(fontsize = fsize)
		       )
rowname_order = rownames(js_df)[row_order(ht_population_js)]

# Get population num and annot
population_num = sub("_.*$", "", clusters_to_plot)
population_annot = sub("^[^_]*_", "", clusters_to_plot)
population_annot_cols = colours_mapping[population_annot]
# Retrieve 
annot_df = data.frame(population_num = population_num, 
		       population_annot = population_annot,
		       row.names = clusters_to_plot)

if(length(heatmaptype)>0){
	if(heatmaptype=='name'){
	annot_df_name = colours_df[colours_df$annotation %in% colnames(js_df),]
	rownames(annot_df_name) = annot_df_name$annotation
	annot_df_name = annot_df_name[rowname_order,]
	annot_df_name$samplename = sub(paste0("^.*?", '_'), "", rownames(annot_df_name))
	# also add TuorLN
	annot_df_name$TuorLN = sub("^(.*?)[0-9].*", "\\1", annot_df_name$samplename)
	annot_df_name$TuorLNColour = colours_mapping[annot_df_name$TuorLN]
	print(annot_df_name)

	row_annot = HeatmapAnnotation(
			   tuln_annot = as.vector(annot_df_name[rowname_order, 'TuorLN'])
			   , ht_annot = as.vector(annot_df_name[rowname_order, 'annotation'])
			   , col = list(ht_annot = setNames(annot_df_name[rowname_order,'colour'], 
							    annot_df_name[rowname_order,'annotation']),
					tuln_annot = setNames(annot_df_name[rowname_order,'TuorLNColour'], 
							      annot_df_name[rowname_order,'TuorLN']))
			   , names = anno_text(annot_df_name[rowname_order, 'samplename'])
			   , show_legend = c(ht_annot = FALSE, tuln_annot = FALSE)
			   , annotation_name_gp = gpar(fontsize = 0)
			   , which = 'row'
			   )
	}
} else {
	row_annot = HeatmapAnnotation(
			   ht_annot = as.vector(annot_df[, 'population_annot'])
			   , clusters = anno_text(annot_df[, 'population_num'], 
					  gp = gpar(fontsize = fsize))
			   , col = list(ht_annot = colours_mapping[annot_df[,'population_annot']])
			   , show_legend = c(ht_annot = FALSE)
			   , annotation_name_gp = gpar(fontsize = 0)
			   , which = 'row'
			   )

}


# Get row annotation
ht_population_js = Heatmap(js_df[rowname_order, rownames(annot_df_name)]
		       , show_row_names = FALSE
		       , row_order = rowname_order
		       , row_names_side = 'right'
		       , show_column_names = FALSE 
		       , cluster_columns = TRUE 
		       , cluster_rows = TRUE 
		       , col = col_mapping
		       , show_heatmap_legend = TRUE 
		       , left_annotation = row_annot
		       , width = 4
		       , row_names_gp = gpar(fontsize = fsize)
		       )

str_to_add = 'SelfMat_noLegend_JS'
filename = get_fname(prefix_filename_for_figures, str_to_add, date)
print(filename)
pdf(filename, width = 5, height = 4)
draw(ht_population_js)
dev.off()
print(paste0('PLOTTED', str_to_add))
