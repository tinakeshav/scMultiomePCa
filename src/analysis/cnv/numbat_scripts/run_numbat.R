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
parser$add_argument('--alleleCounts', type = "character", 
		    required = TRUE, help = "allele counts file")
parser$add_argument('--maxEntropy', type = "double", 
		    required = TRUE, help = "maximum entropy parameter")
parser$add_argument('--nCores', type = "integer", 
		    required = TRUE, help = "number of cores")
parser$add_argument('--multiomeObj', type = "character", 
		    required = TRUE, help = "multiome object for merged patient")
parser$add_argument('--referenceRDSObj', type = "character", 
		    required = TRUE, help = "a reference object as rds format but aggregated")
parser$add_argument('--annotationFile', type = "character", 
		    required = FALSE, help = "annotation csv file")
parser$add_argument('--annotationNameToKeep', type = "character", 
			nargs = '+',
		    required = FALSE, help = "patient name to kee under name col in annotation file")
parser$add_argument('--annotationFileAnnColname', type = "character",  
			nargs = '+', required = FALSE,
		    help = "column name(s) that has the relevant annotation(s)")
parser$add_argument('--annotationFileClusterValue', type = "character", default='wsnn_res.0.4', required = FALSE,
		    help = "cluster to look at under the cluster_col_name of annot file")
args <- parser$parse_args()

print(args)

sample_name = args$sampleName
label = args$labelName
label_dir = args$labelDir
out_dir = args$outDir
max_entropy = args$maxEntropy
ncores = args$nCores
MULTIOME_OBJECT_NAME = args$multiomeObj
referenceRDSObj = args$referenceRDSObj
ann_file = args$annotationFile
name_to_keep = args$annotationNameToKeep
anncol = args$annotationFileAnnColname
annclusterval = args$annotationFileClusterValue



count_mat = c()
df = c()

# Add sample_name to the beginning of 'cell' column in our df
# Such that it matches the barcodes in the multiome object
# and is unique to each sample (avoid duplicated barcodes, although rare)
df[[sample_name]] = fread(glue('{sample_name}/{sample_name}_allele_counts.tsv.gz'), sep = '\t')
df[[sample_name]]$cell = paste0(sample_name, '_', df[[sample_name]]$cell)

sample  = readRDS(MULTIOME_OBJECT_NAME)
normalized_mat_ref = readRDS(referenceRDSObj)

numbat_barcodes = colnames(sample)[colnames(sample) %in% df[[sample_name]]$cell]
count_mat_obs_sample = sample@assays$RNA@counts[,numbat_barcodes]
print('There are', ncol(sample), 'barcodes in the sample')
print('There are', length(unique(df[[sample_name]]$cell)), 'barcodes in the df allele file')
print('There are', ncol(count_mat_obs_sample), 'barcodes in the sample and allele counts file')

out = run_numbat(
		count_mat_obs_sample,
		normalized_mat_ref,
		df[[sample_name]],
		genome = "hg38",
		t = 1e-5,
		ncores = ncores,
		max_entropy = max_entropy,
		multi_allelic = TRUE,
		out_dir = out_dir
)