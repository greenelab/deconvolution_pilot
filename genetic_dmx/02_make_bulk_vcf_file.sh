# Based on https://github.com/lmweber/snp-dmx-cancer, file genotype/genotype_bulk_bcftools/genotype_bulk_HGSOC_bcftools.sh
# This will need to be run once per pool of four samples per type of bulk sequencing, with the pool ID and type of bulk
# sequencing (i.e. dissociated or not, poly-A selected or rRNA depleted) passed in as parameters.

# Options are 12162021 and 01132022
pool=$1

# Options are chunk_ribo, dissociated_ribo, and dissociated_polyA
bulk_type=$2

# Get sample ids based on pool id
if [ $pool == 12162021 ]; then
	sample1=2380
	sample2=2267
	sample3=2283
	sample4=2293
elif [ $pool == 01132022 ]; then
	sample1=2251
	sample2=2428
	sample3=2467
	sample4=2497
else
	echo "Something went wrong"
	exit
fi	

original_location=`pwd`
cd ../../sc-cancer-hgsc/data/bulk_tumors
tumor_location=`pwd`
cd ../index
index_location=`pwd`
cd $original_location

bcftools mpileup -Ou \
	-f $index_location/refdata-gex-GRCh38-2020-A/fasta/genome.fa \
	$tumor_location/$sample1/$bulk_type/STAR/Aligned.sortedByCoord.out.bam \
	$tumor_location/$sample2/$bulk_type/STAR/Aligned.sortedByCoord.out.bam \
	$tumor_location/$sample3/$bulk_type/STAR/Aligned.sortedByCoord.out.bam \
	$tumor_location/$sample4/$bulk_type/STAR/Aligned.sortedByCoord.out.bam | \
bcftools call -mv -Ov \
	-o $tumor_location/pooled_vcf/bcftools_${pool}.vcf

# Create file containing updated sample names
echo "$tumor_location/$sample1/$bulk_type/STAR/Aligned.sortedByCoord.out.bam $sample1" > sample_names_bulk_${pool}_${bulk_type}.txt
echo "$tumor_location/$sample2/$bulk_type/STAR/Aligned.sortedByCoord.out.bam $sample2" >> sample_names_bulk_${pool}_${bulk_type}.txt
echo "$tumor_location/$sample3/$bulk_type/STAR/Aligned.sortedByCoord.out.bam $sample3" >> sample_names_bulk_${pool}_${bulk_type}.txt
echo "$tumor_location/$sample4/$bulk_type/STAR/Aligned.sortedByCoord.out.bam $sample4" >> sample_names_bulk_${pool}_${bulk_type}.txt

# Update sample names in VCF file
bcftools reheader -s sample_names_bulk_${pool}_${bulk_type}.txt $tumor_location/pooled_vcf/bcftools_${pool}_${bulk_type}.vcf > $tumor_location/pooled_vcf/bcftools_${pool}_${bulk_type}_rehead.vcf 
