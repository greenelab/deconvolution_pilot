# Based on https://github.com/lmweber/snp-dmx-cancer, file workflow/scripts/run_cellSNP.sh
# This will need to be run once per pool per bulk sequencing type, with the pool ID and 
# bulk type passed in as parameters.

# Options are 12162021 and 01132022
pool=$1

# Options are chunk_ribo, dissociated_ribo, and dissociated_polyA
bulk_type=$2

original_location=`pwd`
cd ../../sc-cancer-hgsc/data/bulk_tumors/pooled_vcf
vcf_location=`pwd`
cd ../../pooled_tumors/${pool}
pooled_location=`pwd`
cd $original_location

cellsnp-lite \
	-s $pooled_location/bam/pooled.bam \
	-b barcodes_${pool}.txt \
	-O $pooled_location/cellSNP/$bulk_type/ \
	-R $vcf_location/bcftools_${pool}_${bulk_type}_rehead.vcf \
	-p 10 \
	--minMAF=0.1 \
	--minCOUNT=20 \
	--gzip
