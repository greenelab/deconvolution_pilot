# Based on https://github.com/lmweber/snp-dmx-cancer, file workflow/scripts/run_cellSNP.sh
# This will need to be run once per pool, with the pool ID passed in as a parameter.

# Options are 12162021 and 01132022
pool=$1

original_location=`pwd`
cd ../../sc-cancer-hgsc/data/bulk_tumors/pooled_vcf
vcf_location=`pwd`
cd ../../pooled_tumors/${pool}
pooled_location=`pwd`
cd $original_location

cellsnp-lite \
	-s $pooled_location/bam/pooled.bam \
	-b barcodes_${pool}.txt \
	-O $pooled_location/cellSNP \
	-R $vcf_location/bcftools_${pool}_rehead.vcf \
	-p 10 \
	--minMAF=0.1 \
	--minCOUNT=20 \
	--gzip
