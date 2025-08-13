###############################################################################
library(argparse)
library(dplyr)
library(ggplot2)
library(circlize)
library(ComplexHeatmap)
library(data.table)

parser <- ArgumentParser(description='Plotting entropy in ME versus cell cycle')

parser$add_argument('--metadata_filenames',
                    type = "character",
                    nargs = '+',
                    required = TRUE,
                    help = "metadata filename to read and \\
                    generate column medians for; has to \\
                    include all samples in the same order as the sample_entropy_filenames")

parser$add_argument('--sample_names',
                    type = "character",
                    nargs = '+',
                    required = TRUE,
                    help = "sample names to add as identified and grab colours from")

parser$add_argument('--isolate_ME', 
                    type = "logical", 
                    required = FALSE, 
                    default = TRUE,
                    help = "Set to FALSE to note isolate ME")

parser$add_argument('--max_ent_on_nsamples', 
                    type = "logical", 
                    required = FALSE, 
                    default = TRUE,
                    help = "If TRUE then it computes max entropy based on metadata sample column name")

parser$add_argument('--date', 
                    type = "character",
                    required = TRUE, 
                    help = "date for filenames that get saved")

args <- parser$parse_args()
print(args)

metadata_filenames = args$metadata_filenames
sample_names = args$sample_names
sample_entropy_filenames = args$sample_entropy_filenames
isolate_ME = args$isolate_ME
max_ent_on_nsamples = args$max_ent_on_nsamples
date = args$date

###############################################################################
# Functions
###############################################################################
entropy <- function(x) {
  prob <- table(x) / length(x)
  -sum(prob * log2(prob))
}

entropy_on_metadata <- function(metadata_df){
    x_entropy_col = 'cluster_annot_short'
    y_entropy_cols = 'name'
    entropy_df = data.table(metadata_df)[,lapply(.SD, entropy), by = x_entropy_col, .SDcols = y_entropy_cols]
    setnames(entropy_df, c(x_entropy_col, y_entropy_cols), c('population', 'name'))
    return(as.data.frame(entropy_df))
}

compute_max_entropy <- function(num_categories) {
    return(log(num_categories) / log(2))  # Return entropy in bits
}

normalize_entropy_values <- function(values, max_value) {
    min_value <- 0  # Set minimum to 0
    normalized <- (values - min_value) / (max_value - min_value)  # Normalize
    return(normalized)
}

get_identifier_colours <- function(){
    row_annot_cols_filename = 'data/colours_23032023_with_samples_averaged_vals_03032023_CNVrep.csv'
    row_annot_col_df = read.csv(row_annot_cols_filename)
    row_annot_list = setNames(row_annot_col_df$colour,
                            row_annot_col_df$annotation)

    # can skip this, its for other purposes
    # Get population num and annot
    #population_num = sub("_.*$", "", rowname_order)
    #population_annot = sub("^[^_]*_", "", rowname_order)
    #population_annot_cols = row_annot_list[population_annot]
    #annot_df = data.frame(population_num = population_num,
    #                    population_annot = population_annot,
    #                    row.names = rowname_order)
    return(row_annot_list)
}


get_population_annotations_and_colours <- function(entropy_df){
    row_annot_cols_filename = 'data/colours_23032023_with_samples_averaged_vals_03032023_CNVrep.csv'
    row_annot_col_df = read.csv(row_annot_cols_filename)
    row_annot_list = setNames(row_annot_col_df$colour,
                            row_annot_col_df$annotation)
    population_cancer_list = setNames(row_annot_col_df$cancer,
                                      row_annot_col_df$annotation)

    # can skip this, its for other purposes
    # Get population num and annot
    entropy_df$population_num = sub("_.*$", "", entropy_df$population)
    entropy_df$population_annot = sub("^[^_]*_", "", entropy_df$population)
    entropy_df$population_annot_cols = row_annot_list[entropy_df$population_annot]

    # add population cancer versus tme
    entropy_df$is_cancer = population_cancer_list[entropy_df$population_annot]

    #annot_df = data.frame(population_num = population_num,
    #                      population_annot = population_annot,
    #                     row.names = rowname_order)
    #return(row_annot_list)
    return(entropy_df)
}

############################################################
# Get colours
colours = get_identifier_colours()
print(head(colours))
# Initialize lists to store medians and maximum entropies
max_entropies = list()
entropy_data = list()
# Loop through both metadata and sample entropy filenames

for (i in seq_along(metadata_filenames)) {
    # Read metadata file and calculate column (cell cycle module score) medians
    metadata_data <- read.csv(metadata_filenames[i])  
    
    median_data = metadata_data %>% 
                  group_by(cluster_annot_short) %>% 
                  summarise(median_cellcycle = median(G2_M_markers1)) %>%
                  as.data.frame() 
    rownames(median_data) = median_data$cluster_annot_short

    # Read sample entropy file and calculate maximum entropy
    #entropy_data[[i]] <- read.csv(sample_entropy_filenames[i])
    entropy_data[[i]] = entropy_on_metadata(metadata_data)
    rownames(entropy_data[[i]]) = entropy_data[[i]]$population

    # Calculate maximum entropy based on number of categories
    max_entropies[[i]] <- compute_max_entropy(length(unique(metadata_data$name))) 

    # Normalize entropies
    normalized_entropies = normalize_entropy_values(values = entropy_data[[i]][['name']], 
                                                        max_value = max_entropies[[i]])
    entropy_data[[i]]$normalized_entropy = normalized_entropies

    # Get name and colour for entropy data
    entropy_data[[i]]$identifier = rep(sample_names[i], 
                                       nrow(entropy_data[[i]]))
    entropy_data[[i]]$identifier_colour = rep(colours[[sample_names[i]]], 
                                              nrow(entropy_data[[i]]))
    # Add medians per column
    entropy_data[[i]]$G2_M_markers = median_data[rownames(entropy_data[[i]]), 'median_cellcycle']

    # Add population annotations
    entropy_data[[i]] = get_population_annotations_and_colours(entropy_data[[i]])

    
}

# Bind entropy data to plot
entropy_data_full = do.call(rbind, entropy_data)
print(entropy_data_full)
filename_entropy_data_full = paste0('data/allPatsTuONLY_EntropyNormMat_', date, '.csv')
write.csv(entropy_data_full, filename_entropy_data_full)


####################################################################
# Plot annotations / entropy / cell cycle heatmaps 
####################################################################
patient_plots = list()
for (i in seq_along(max_entropies)) {
    # patient we will be plotting
    entropy_df = entropy_data[[i]] %>% arrange(is_cancer, desc(name))

    # entropy colour per patient identifier
    entropy_col = colorRamp2(c(0,max_entropies[[i]]),
                             c('#F4E8E6', "#074517"))
    
    # manually obtained to cover entire range
    g2m_col = colorRamp2(c(-0.1, 0, 0.5),
                         c("#7B68CC", "#F4E8E6", "#F7A35C"))

    # annotation colours
    annot_col = setNames(entropy_df$population_annot_cols,
                         entropy_df$population_annot)

    row_a = HeatmapAnnotation(
            clusters = anno_text(entropy_df[,'population_num'], rot = 0, just = 'centre') 
            , ht_annot = as.vector(entropy_df[,'population_annot'])
            , ent_annot = entropy_df[, 'name']
            , cellcycle_annot = entropy_df[, 'G2_M_markers']
            , col = list(ht_annot = annot_col 
                            , ent_annot = entropy_col
                            , cellcycle_annot = g2m_col)
            , show_legend = c(ht_annot = FALSE)
            , annotation_name_gp = gpar(fontsize = 0)
            , which = 'column'
            )
    patient_plots[[i]] = row_a  
    filename_combined = paste0('figures/', 
                               sample_names[i],
                               '_combined_heatmap_maxEntropy_CellCycle_', 
                               date, '.pdf')
    pdf(filename_combined, height = 1, width = 8)
    draw(row_a)
    dev.off()
    print(filename_combined)
}

# save cell cycle legend independently 
filename_legend = paste0('figures/', 
                        'allPats_legendCellCycle__', 
                        date, '.pdf')
pdf(filename_legend, height = 3, width = 4)
lgd = Legend(col_fun = g2m_col, title = "Cell Cycle Score")
draw(lgd)
dev.off()

####################################################################
# Plot correlation plot  
####################################################################
if (isolate_ME) {
    entropy_data_me <- entropy_data_full[grepl('_Malignant Epithelial', entropy_data_full$population), ]
}

p = ggplot(entropy_data_me, 
           aes(x = normalized_entropy, y = G2_M_markers, color = G2_M_markers)) +
    geom_point() +
    labs(x = "Normalized entropy",
         y = "Cell cycle score in malignant populations",
         title = "Scatter Plot of Normalized Maximum Entropy vs Median Correlation") +
    theme_classic() +
    scale_color_gradientn(colors = c("#7B68CC", "#F4E8E6", "#F7A35C"))

filename = paste0('figures/allPats_ScatterPlot_maxEntropy_CellCycle_', date, '.pdf')
ggsave(filename, plot = p, height = 3, width = 5)
