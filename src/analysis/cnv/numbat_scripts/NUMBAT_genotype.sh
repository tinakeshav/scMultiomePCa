#!/bin/bash


output_folder="$1"
label="$2"
sample_name="$3"
shift 3
vcf_base=("$@")

mkdir -p $output_folder

Rscript multiome_scripts/NUMBAT_genotype.R \
	--output_folder $output_folder \
	--sample_name $sample_name \
	--label $label \
	--vcf_base "${vcf_base[@]}"

# sample_name is never used in the numbat genotyping function.

# Check if the uncompressed chr1 VCF exists but the compressed version and index don't
if [ -f "${output_folder}/${label}_chr1.vcf" ] && [ ! -f "${output_folder}/${label}_chr1.vcf.gz" ] && [ ! -f "${output_folder}/${label}_chr1.vcf.gz.tbi" ]; then
    input_dir=$(ls ${output_folder}/*chr*.vcf)

    for vcffile in $input_dir
    do
        bgzip $vcffile
        tabix $vcffile'.gz'
    done
    echo 'Finished running bgzip and tabix'
else
    echo "Skipping bgzip and tabix: Either ${output_folder}/${label}_chr1.vcf doesn't exist or compressed files already exist."
fi
