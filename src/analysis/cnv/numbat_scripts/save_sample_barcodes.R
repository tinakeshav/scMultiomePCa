# save 
# save the tumor barcodes for each sample
# for future runs on patient-level analysis because
# not all cells from the patient fit in memory :/ 

library(data.table)
library(dplyr)
library(argparse)


# Get arguments
parser <- ArgumentParser(description='Save sample barcodes, subsampled to tumor barcodes / in populations of interest')

parser$add_argument('--clone_post_2', type = "character", 
		    required = TRUE, help = "clone_post_2.tsv file")

parser$add_argument('--annotation_file', type = "character", 
		    required = TRUE, help = "annotation file")

parser$add_argument('--annotation_colname', type = "character", 
            default = "cluster_annot_short",
		    required = FALSE, help = "annotation column name")

parser$add_argument('--annotation_cluster_value', type = "character", 
            default = "_Malignant Epithelial",
		    required = FALSE, help = "annotation cluster value")

parser$add_argument('--out_dir', type = "character", 
        default = "./",
	    required = FALSE, help = "output directory")

args <- parser$parse_args()

clone_post_2_path <- args$clone_post_2
annotation_file_path <- args$annotation_file
annotation_colname <- args$annotation_colname
annotation_cluster_value <- args$annotation_cluster_value
out_dir <- args$out_dir

# read numbat clonal calls / tumor vs normal comparmetns 
clone_dt <- fread(clone_post_2_path)
tumor_cells <- clone_dt %>%
  filter(compartment == "tumor") %>%
  pull(cell) %>%
  unique()

# read metadata with annotations
ann_dt <- fread(annotation_file_path)
ann_dt[, cell_barcode := V1]
ann_dt$cell = ann_dt$barcode

if (!(annotation_colname %in% colnames(ann_dt))) {
  stop(paste0("Annotation column ", annotation_colname, " not found in annotation file"))
}

# filter df to only contain annotation rows where 
# annotation_colname contains the exact annotation_cluster_value string
mask <- !is.na(ann_dt[[annotation_colname]]) & grepl(annotation_cluster_value, ann_dt[[annotation_colname]], fixed = TRUE)
ann_filtered <- ann_dt[mask]
ann_barcodes <- ann_filtered$cell_barcode

# intersect tumor cells from clone_post_2 with annotation barcodes
selected_barcodes <- intersect(ann_barcodes, tumor_cells)
# If length of selected_barcodes is less than tumorBarcodesSubsample, 
# keep everything. If not, subsample to tumorBarcodesSubsample
set.seed(22)
if (length(selected_barcodes) < tumorBarcodesSubsample) {
  selected_barcodes <- selected_barcodes
} else {
  selected_barcodes <- sample(selected_barcodes, tumorBarcodesSubsample)
}

# write barcodes to txt/tsv file
outfile <- file.path(out_dir, "tumor_barcodes.txt")
barcodes_to_write = data.table(cell = selected_barcodes)
fwrite(barcodes_to_write, outfile, col.names = TRUE)
message(sprintf("Saved %d tumor barcodes to %s", length(selected_barcodes), outfile))