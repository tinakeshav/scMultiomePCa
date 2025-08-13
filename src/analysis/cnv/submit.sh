#!/bin/bash
#SBATCH --mem 2GB
#SBATCH  -c 1
#SBATCH  -t 24:00:00
#SBATCH -p all
#SBATCH -A cbmp
#SBATCH --output=R-%x.%j.out
#SBATCH --error=R-%x.%j.err

## inside scvi
snakemake --snakefile Snakefile_NUMBAT -c 1 -j 60  --cluster "sbatch {cluster.params}" -u cluster_config.yaml --latency-wait 300 --config samples_NUMBAT.csv --dryrun
snakemake --snakefile Snakefile_PATIENT_multisample -c 1 -j 60  --cluster "sbatch {cluster.params}" -u cluster_config.yaml --latency-wait 300 --config samples_file=samples_PATIENT_multisamples_Patient3.csv --dryrun
