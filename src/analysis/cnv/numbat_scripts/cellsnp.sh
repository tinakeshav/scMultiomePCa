#!/bin/bash

### ensure that the environment has cellsnp-lite
input_bam=$1
input_barcodes=$2
output_directory=$3
snpvcf=$4
cores=$5
umitag=$6

echo $input_bam
echo $input_barcodes
echo $output_directory
echo $snpvcf
echo $cores

mkdir -p $output_directory

rel_path=$(realpath --relative-to $(pwd) $output_directory)
echo $rel_path

cellsnp-lite -s $input_bam -b $input_barcodes -O $rel_path -R $snpvcf -p $cores --minMAF 0 --minCOUNT 2 --UMItag $umitag --cellTAG CB 
