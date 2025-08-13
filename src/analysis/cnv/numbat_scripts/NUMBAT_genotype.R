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

print(getRversion())
parser <- ArgumentParser(description='Run SNP pileup and phasing with 1000G')

parser$add_argument('--output_folder', 
                    type = "character", 
                    required = TRUE, 
                    help = "Individual output folder with phasing folder already inside")

parser$add_argument('--sample_name', 
                    type = "character", 
                    required = FALSE, 
                    help = "sample name to pass to numbat::genotype")

parser$add_argument('--vcf_base', 
                    type = "character", 
		    nargs = '+',
                    required = TRUE, 
                    help = "vcf files")

parser$add_argument('--label', 
                    type = "character", 
                    required = FALSE, 
                    help = "Label folder to save chromosome-based vcf files. For multi-input phasing")

args <- parser$parse_args()
print(args)

sample_name = args$sample_name
label = args$label
vcf_base = args$vcf_base
output_folder = args$output_folder

## VCF creation
cat('Creating VCF\n')
#vcf = vcfR::read.vcfR(vcf_base, verbose = F)

vcfs = lapply(vcf_base, function(vcf_file) {
	print(vcf_file)
	if (file.exists(vcf_file)) {
		if (file.size(vcf_file) != 0) {
			vcf = vcfR::read.vcfR(vcf_file, verbose = F)
			# Remove chr prefix if present -- from numbat pileup and phasing script
			vcf@fix[,1] <- gsub("chr", "", vcf@fix[,1])
			message('Read pileup VCF')
			return(vcf)
		} else {
			stop('Pileup VCF is empty')
		}
		} else {
			stop('Pileup VCF not found')
		}
})


cat('Running Numbat Genotype on VCF\n')


# Run numbat genotype with the processed VCF
# numbat:::genotype's second argument (samples) is never actually used in the function body
# can be anything...?
numbat:::genotype(label, sample_name, vcfs, output_folder, chr_prefix = TRUE)
