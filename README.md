# scMultiome Analyses
Data required for directly reproducing the figures are uploaded on EGA and GEO.
The raw data can be fed to the Snakemake files for various steps in the analyses ```src/workflow/Snake*``` 
- ```Snakefile_MULTIOME_preAnnotation``` : Generate multiome files from raw 10X outputs (on GEO) and runs WNN analyses
- ```Snakefile_CellBender``` : Runs CellBender on 10X outputs with sample-specific parameters. Outputs can be used to re-run the above Snakemake file as well.


Analyses and plotting functions 
- ```entropy``` 
- ```TF chromvar enrichment```
- ```Pythonic inter-patient```
- ```Intra-patient correlations```
- ```CNV analysis and plots```
- ```comparisons of cellbender```


All raw FASTQs are accessed through EGA study: EGAS50000000524.
All pre-processed and processed multiome and metadata files used in all scripts addressed here can be downloaded from GEO: GSE271675.
