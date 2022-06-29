# Based on https://github.com/lmweber/snp-dmx-cancer, file genotype/align_index_bulk_STAR/align_index_bulk_HGSOC_17667X1.sh
# This will need to be run once for each sample and for each type of bulk sequencing, with the sample ID and the bulk type
# (i.e. dissociated or not, rRNA depleted or poly-A selected) passed in as parameters.

# Options are 2251, 2267, 2283, 2293, 2380, 2428, 2467, 2497
sample=$1
# Options are chunk_ribo, dissociated_ribo, and dissociated_polyA 
bulk_type=$2

original_location=`pwd`
cd ../../sc-cancer-hgsc/data/index
index_location=`pwd`
cd ../bulk_tumors/$sample/$bulk_type
tumor_location=`pwd`

cd $original_location

# STAR only accepts one set of files, so we'll merge across lanes
cat $tumor_location/${sample}_L001_R1_001.fastq.gz $tumor_location/${sample}_L002_R1_001.fastq.gz > $tumor_location/${sample}_R1_merged.fastq.gz
cat $tumor_location/${sample}_L001_R2_001.fastq.gz $tumor_location/${sample}_L002_R2_001.fastq.gz > $tumor_location/${sample}_R2_merged.fastq.gz

STAR \
	--genomeDir $index_location/STAR_index \
	--runThreadN 6 \
	--readFilesIn $tumor_location/${sample}_R1_merged.fastq.gz  $tumor_location/${sample}_R2_merged.fastq.gz \
	--outFileNamePrefix $tumor_location/STAR/ \
	--readFilesCommand gunzip -c \
	--outSAMtype BAM SortedByCoordinate \
	--quantMode GeneCounts

samtools index $tumor_location/STAR/Aligned.sortedByCoord.out.bam
