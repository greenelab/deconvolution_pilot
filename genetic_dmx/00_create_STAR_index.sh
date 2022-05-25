# Based on https://github.com/lmweber/snp-dmx-cancer, file genotype/align_index_bulk_STAR/create_STAR_index.sh
# This should only need to be run once, so it will check if the file exists first.

current_location=`pwd`
cd ../../sc-cancer-hgsc/data/index
index_location=`pwd`
cd $current_location

# Check if file already exists
if [ -f "$index_location/STAR_index/SAindex" ]; then
	echo "STAR index exists, skipping."
else
	STAR \
		--runMode genomeGenerate \
		--runThreadN 6 \
		--genomeDir $index_location/STAR_index \
		--genomeFastaFiles $index_location/refdata-gex-GRCh38-2020-A/fasta/genome.fa \
		--sjdbGTFfile $index_location/refdata-gex-GRCh38-2020-A/genes/genes.gtf
fi
