# copy of run-numbat_afterPhasing_AnyPat_Clonal.R but supplies a reference as rds file
library(numbat)
library(dplyr)
library(data.table)
library(glue)
library(stringr)
library(Matrix)
library(magrittr)
library(parallel)
library(Seurat)
library(Signac)
library(argparse)

parser <- ArgumentParser(description='Patient Clonal samples')

parser$add_argument('--sampleName', type = "character", 
		    required = TRUE, help = "sample name", nargs = '+')

parser$add_argument('--alleleCounts', type = "character", nargs = '+',
		    required = TRUE, help = "allele counts file")

parser$add_argument('--tumorBarcodes', type = "character", 
		    required = TRUE, help = "tumor barcodes file", nargs = '+')

parser$add_argument('--tumorBarcodesSubsample', type = "double", 
		    required = TRUE, help = "number of tumor barcodes to subsample to from each sample")

parser$add_argument('--labelName', type = "character", 
		    required = FALSE, help = "label name")

parser$add_argument('--labelDir', type = "character", 
		    required = FALSE, help = "label dir with all the sample allele counts in it")

parser$add_argument('--outDir', type = "character", 
		    required = TRUE, help = "out directory to save results to")

parser$add_argument('--maxEntropy', type = "double", 
		    required = TRUE, help = "maximum entropy parameter")

parser$add_argument('--nCores', type = "integer", 
		    required = TRUE, help = "number of cores")

parser$add_argument('--multiomeObj', type = "character", 
		    required = TRUE, help = "multiome object for merged patient")

parser$add_argument('--referenceRDSObj', type = "character", 
		    required = TRUE, help = "a reference object as rds format but aggregated")

args <- parser$parse_args()

print(args)

sample_names = args$sampleName
label = args$labelName
label_dir = args$labelDir
out_dir = args$outDir
max_entropy = args$maxEntropy
ncores = args$nCores
MULTIOME_OBJECT_NAME = args$multiomeObj
referenceRDSObj = args$referenceRDSObj
tumorBarcodes = args$tumorBarcodes
tumorBarcodesSubsample = args$tumorBarcodesSubsample


sample  = readRDS(MULTIOME_OBJECT_NAME)
reference = readRDS(referenceRDSObj)

count_mat = c()
df = c()
subsampled_tumor_barcodes = rbindlist(lapply(tumorBarcodes, fread))


# Add sample_name to the beginning of 'cell' column in our df
# Such that it matches the barcodes in the multiome object
# and is unique to each sample (avoid duplicated barcodes, although rare)
for (sample_name in length(sample_names)) {
    df[[sample_name]] = fread(glue('{label_dir}/{sample_name}_allele_counts.tsv.gz'), sep = '\t')
    barcodes[[sample_name]] = df[[sample_name]]$cell
	df[[sample_name]]$cell = paste0(sample_name, '_', df[[sample_name]]$cell)

}

normalized_mat_ref = readRDS(referenceRDSObj)
count_mat_obs_sample = sample@assays$RNA@counts[,subsampled_tumor_barcodes]

out = run_numbat(
		count_mat_obs_sample,
		normalized_mat_ref,
		df[sample_names] %>% bind_rows %>% dplyr::filter(cell %in% subsampled_tumor_barcodes),
		genome = "hg38",
		t = 1e-5,
		ncores = ncores,
		max_entropy = max_entropy,
		multi_allelic = TRUE,
		out_dir = out_dir
)