# Based on https://github.com/lmweber/snp-dmx-cancer, file workflow/scripts/run_vireo.sh
# This will need to be run once per pool, with the pool ID passed in as a parameter.

# Options are 12162021 and 01132022
pool=$1

original_location=`pwd`
cd ../../sc-cancer-hgsc/data/pooled_tumors/${pool}
pooled_location=`pwd`
cd $original_location

vireo -c $pooled_location/cellSNP -N 4 -o $pooled_location/vireo --randSeed=123
