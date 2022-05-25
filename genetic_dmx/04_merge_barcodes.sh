# cellSNP requires a list of all barcodes used; we'll borrow this from the cellranger multi output.
# This will need to be run once per pool of four samples, with the pool ID manually entered on line 4.

# Options are 12162021 and 01132022
pool=01132022

original_location=`pwd`
cd ../../sc-cancer-hgsc/data/pooled_tumors/${pool}
tumor_location=`pwd`

cd $original_location

cut $tumor_location/Cellranger/outs/multi/multiplexing_analysis/assignment_confidence_table.csv -d',' -f6 | sed '1d' > barcodes_${pool}.txt
