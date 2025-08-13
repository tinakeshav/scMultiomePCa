#!/bin/bash

# The directory containing the fragment files to call peaks on
frag_dir=$1
# Where macs results are saved in 
macs_dir=$2
# The suffix to fragment files (since prefix is 'unknown', cluster+annot)
bed_file_suffix=$3
# The output text file to save results in 
output_txtfile=$4

files=$(ls $frag_dir/*$bed_file_suffix)

mkdir -p $macs_dir

for file in $files
	do
		samplename="${file##*/}"
		samplename="${samplename%.bed}"
		echo $samplename
		macs3 callpeak -t $file -n $samplename --outdir $macs_dir -f BED -q 0.01 -g 2.7e9 --keep-dup all --nomodel --nolambda --shift -75 --extsize 150 --bdg --call-summits
		if [ $? -ne 0 ]; then
			echo "macs3 callpeak failed for $file"
			exit 1
		fi
		# Check for expected output file (narrowPeak)
		peak_file="$macs_dir/${samplename}_peaks.narrowPeak"
		if [ ! -f "$peak_file" ]; then
			echo "Expected output $peak_file not found for $file"
			exit 1
		fi
	done
# NOTE: never tested on failed files. 
# It always worked. So not sure if it'll proceed THERE IS NOTHING TAKING CARE OF FAILURE HERE. 
# Only create output_txtfile if all macs3 calls succeeded

touch $output_txtfile
