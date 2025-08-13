# script to call intra-patient peaks and compute fragments counts
# Purpose and associated LBN page: 
# FIGURE 1; SCRIPT : ATAC non-coding and promoter region correlations

###################################################################
# Load Argparse
###################################################################
suppressPackageStartupMessages(library("argparse"))

parser = ArgumentParser()

parser$add_argument("--multiomeobjname", 
		    type = 'character',
		    required = TRUE,
		    help="multiome object to load. ")

parser$add_argument("--metadatatoattach", 
		    type = 'character',
		    required = TRUE,
		    help="metadata file to attach")

parser$add_argument("--colname_for_frags", 
		    type = 'character',
		    required = TRUE,
		    help="column for macs grouping")

parser$add_argument("--date", 
		    type="character",
		    required = TRUE,
		    help="date script is running ; \\
		    will be used in saving object names")

parser$add_argument("--frag_dir", 
		    type="character",
		    help="fragment directory to save files in")

parser$add_argument("--frag_suffix", 
		    type="character",
		    help="suffix before .bed")

parser$add_argument("--output_txtfile", 
		    type="character",
		    required = TRUE,
		    help="output txtfile name to write ; \\
		    after asserting whether all files are present")

args <- parser$parse_args()

print(args)

multiome_obj_name =  args$multiomeobjname
metadatatoattach = args$metadatatoattach
colname_for_frags = args$colname_for_frags
date = args$date
frag_dir = args$frag_dir
frag_suffix = args$frag_suffix
output_txtfile = args$output_txtfile

################################################################################
# Load the rest of the required libraries 
# Plus function for creating output file
################################################################################
suppressMessages(library(Seurat))
suppressMessages(library(Signac))
suppressMessages(library(EnsDb.Hsapiens.v86))
suppressMessages(library(GenomicRanges))


check_frags <- function(directory_path, 
			unique_col_vals, 
			frag_suff,
			output_txtfile) {

	expected_files = paste0(gsub(' ', '_', unique_col_vals), 
				frag_suff)
	expected_files = paste0(expected_files, '.bed')
	expected_file_paths = file.path(directory_path, expected_files)
	print(expected_file_paths)

	# List files in the specified directory
	files_in_directory <- list.files(directory_path, 
					 full.names = TRUE)

	# Check if each expected file exists in the directory
	exp_exists  <-  sapply(expected_file_paths, function(file) {
			file %in% files_in_directory
			})
	print(exp_exists)

	if(all(exp_exists)) {
		file.create(output_txtfile)
		print('attempted to create output fin file')
		return(TRUE)
	}
	return(FALSE)
}

################################################################################
# Load the rest of the required libraries 
# Plus function for creating output file
################################################################################
# Before loading, let's check if output file already exists (indicates the 
# script has been run once before and files are there)

print('checking if necessary to regenerate fragment files')
metadata = read.csv(metadatatoattach, 
		    header = TRUE, 
		    row.names = 1)

res = check_frags(frag_dir, 
		  unique(metadata[, colname_for_frags]), 
		  frag_suffix,
		  output_txtfile)

if(!(res)){
	print('rerunning script')
	object = readRDS(multiome_obj_name)

	# attach metadata
	object = AddMetaData(object, 
			     metadata)

	features = SplitFragments(
			     object,
			     group.by = colname_for_frags,
			     file.suffix = frag_suffix,
			     outdir = frag_dir,
			     )

	# function to check if all bed files are there
	# if pass then create output


	# check if all frag files are present, and if so 
	# create the 'output' txt file
	check_frags(frag_dir, 
		    unique(object@meta.data[, colname_for_frags]), 
		    output_txtfile)
}
