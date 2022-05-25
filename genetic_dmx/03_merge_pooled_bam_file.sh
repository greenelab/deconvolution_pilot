# When we run cellranger multi, it splits up the bam files: one per sample with aligned reads from 
# cells it has assigned to that sample, and one with unassigned alignments. We merge those into one
# file for cellSNP/vireo.
# This will need to be run once per pool of four samples, with the pool ID passed in as a parameter.

# Options are 12162021 and 01132022
pool=$1

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
cd ../../sc-cancer-hgsc/data/pooled_tumors/${pool}/Cellranger/outs
tumor_location=`pwd`
cd ../..
out_location=`pwd`
cd $original_location

samtools merge $out_location/bam/pooled.bam \
       $tumor_location/multi/count/unassigned_alignments.bam \
       $tumor_location/per_sample_outs/${sample1}/count/sample_alignments.bam \
       $tumor_location/per_sample_outs/${sample2}/count/sample_alignments.bam \
       $tumor_location/per_sample_outs/${sample3}/count/sample_alignments.bam \
       $tumor_location/per_sample_outs/${sample4}/count/sample_alignments.bam

samtools index $out_location/bam/pooled.bam
