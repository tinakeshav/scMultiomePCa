#!/usr/bin/env Rscript
suppressPackageStartupMessages({
	library(argparse)
	library(glue)
	library(stringr)
	library(data.table)
	library(dplyr)
	library(vcfR)
	library(Matrix)
	library(numbat)
})


parser <- ArgumentParser(description='Run SNP pileup and phasing with 1000G')
parser$add_argument('--phased_vcfs', type = "character", 
		    required = TRUE, nargs = '+', 
		    help = "All phased VCF files for this run")
parser$add_argument('--pileup_dirs', type = "character", 
		    required = TRUE, nargs = '+', 
		    help = "All pileup directories for this run.")
parser$add_argument('--allele_counts_names', type = "character", 
		    required = TRUE, nargs = '+', 
		    help = "Names for allele counts. eg. {sample}_{assay}")
parser$add_argument('--label', type = "character", 
		    required = TRUE, 
		    help = "Label name used for phasing/genotyping")
parser$add_argument('--outdir', type = "character", 
		    required = TRUE, 
		    help = "Output directory with pileup (sample/assay)-specific allele counts.")
parser$add_argument('--genome', type = "character", 
		    required = TRUE, help = "hg38 or hg19")
parser$add_argument('--gmap', type = "character", 
		    required = TRUE, 
		    help = "Path to genetic map provided by Eagle2 (e.g. Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz)")
args <- parser$parse_args()

print(args)

phased_vcfs = args$phased_vcfs
pileup_dirs = args$pileup_dirs
allele_counts_names = args$allele_counts_names
outdir = args$outdir
label = args$label
gmap = args$gmap
genome = args$genome

## Generate allele count dataframe
cat('Generating allele count dataframes\n')

if (genome == 'hg19') {
    gtf = gtf_hg19
} else {
    gtf = gtf_hg38
}

genetic_map = fread(gmap) %>% 
    setNames(c('CHROM', 'POS', 'rate', 'cM')) %>%
    group_by(CHROM) %>%
    mutate(
        start = POS,
        end = c(POS[2:length(POS)], POS[length(POS)])
    ) %>%
    ungroup()

print(genetic_map)

# read in phased VCF
vcf_phased = lapply(phased_vcfs, function(vcf_file) {
	message('Reading', vcf_file)
	if (file.exists(vcf_file)) {
		fread(vcf_file, skip = '#CHROM') %>%
		rename(CHROM = `#CHROM`) %>%   
		mutate(CHROM = str_remove(CHROM, 'chr'))
	} else {
		stop('Phased VCF not found')
	}
	}) %>%
	Reduce(rbind, .) %>%
	mutate(CHROM = factor(CHROM, unique(CHROM)))

for (i in 1:length(pileup_dirs)){
	allele_counts_name=allele_counts_names[[i]]
	pu_dir=pileup_dirs[[i]]

	# pileup VCF
	message('Reading from ', pu_dir)
	vcf_pu = fread(glue('{pu_dir}/cellSNP.base.vcf'), skip = '#CHROM') %>% 
		rename(CHROM = `#CHROM`) %>%
		mutate(CHROM = str_remove(CHROM, 'chr'))

	# count matrices
	AD = readMM(glue('{pu_dir}/cellSNP.tag.AD.mtx'))
	DP = readMM(glue('{pu_dir}/cellSNP.tag.DP.mtx'))

	cell_barcodes = fread(glue('{pu_dir}/cellSNP.samples.tsv'), header = F) %>% pull(V1)

	# setting sample as label is normal here. Just bad naming of 'sample' <> 'label'; in case of multi-sample
	# it should be the label however.
	# it's referring to the column with phasing information 
	df = numbat:::preprocess_allele(
		sample = label,
		vcf_pu = vcf_pu,
		vcf_phased = vcf_phased,
		AD = AD,
		DP = DP,
		barcodes = cell_barcodes,
		gtf = gtf,
		gmap = genetic_map
		) %>%
		filter(GT %in% c('1|0', '0|1'))
	message('Saving allele counts for sample ', allele_counts_name)

	fwrite(df, glue('{outdir}/{allele_counts_name}_allele_counts.tsv.gz'), sep = '\t')

}
