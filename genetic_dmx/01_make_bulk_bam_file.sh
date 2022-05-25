# Based on https://github.com/lmweber/snp-dmx-cancer, file genotype/align_index_bulk_STAR/align_index_bulk_HGSOC_17667X1.sh
# This will need to be run once for each sample, with the sample ID manually entered in on line 5.

# Options are 2251, 2267, 2283, 2293, 2380, 2428, 2467, 2497
sample=2267

original_location=`pwd`
cd ../../sc-cancer-hgsc/data/index
index_location=`pwd`
cd ../bulk_tumors/$sample/chunk_ribo
tumor_location=`pwd`

# Each sample has a unique run id, e.g. S62 in ${sample}_${SID}_L001_R1_001.fastq.gz
# We'll pull it out here to consistently name the merged files
SIDtmp=`ls *_L001_R1_001.fastq.gz` 
SID="$(echo $SIDtmp | cut -d'_' -f2)"

cd $original_location

# STAR only accepts one set of files, so we'll merge across lanes
cat $tumor_location/${sample}_${SID}_L001_R1_001.fastq.gz $tumor_location/${sample}_${SID}_L002_R1_001.fastq.gz > $tumor_location/${sample}_${SID}_R1_merged.fastq.gz
cat $tumor_location/${sample}_${SID}_L001_R2_001.fastq.gz $tumor_location/${sample}_${SID}_L002_R2_001.fastq.gz > $tumor_location/${sample}_${SID}_R2_merged.fastq.gz

STAR \
	--genomeDir $index_location/STAR_index \
	--runThreadN 6 \
	--readFilesIn $tumor_location/${sample}_${SID}_R1_merged.fastq.gz  $tumor_location/${sample}_${SID}_R2_merged.fastq.gz \
	--outFileNamePrefix $tumor_location/STAR/ \
	--readFilesCommand gunzip -c \
	--outSAMtype BAM SortedByCoordinate

samtools index $tumor_location/STAR/Aligned.sortedByCoord.out.bam
