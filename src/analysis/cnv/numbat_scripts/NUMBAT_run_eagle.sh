#comes from RScript 
#make-shscripts-from-pileup_and_phase.R
eagle=$1
numThreads=$2
vcfChrDir=$3
geneticMapFile=$4
vcfRef=$5
#outSuffix=$5 # goes inside vcf Chr Dir

echo $eagle
echo $numThreads
echo $vcfChrDir
echo $geneticMapFile
echo $vcfRef
#echo $outSuffix # goes inside vcf Chr Dir

for vcfTarget in $(ls $vcfChrDir/*_chr*.vcf.gz)
do
	outName=$(echo $vcfTarget | sed "s/vcf.gz.*/phased/g")
	echo $outName

	# making the appropriate genotype bcf filename
	chrpaneldir_suff=$(echo $outName | awk '{print $NF}' FS=_) # returns chr10.phased
	#echo $chrpaneldir_suff
	chrpaneldir_suff=${chrpaneldir_suff/phased/genotypes.bcf}
	#echo $chrpaneldir_suff
	chrpaneldir=$vcfRef'/'$chrpaneldir_suff
	echo $chrpaneldir

	bash -c "$eagle'eagle' --numThreads $numThreads --vcfTarget $vcfTarget --geneticMapFile $geneticMapFile --outPrefix $outName --vcfRef $chrpaneldir"
done