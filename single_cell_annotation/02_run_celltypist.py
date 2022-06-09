# See https://www.celltypist.org/ for more info

import sys
import celltypist

# This can be run on an individual sample (valid ids 2251, 2267, 2283, 2293,
# 2380, 2428, 2467, 2497), an individual pooled run (valid ids 12162021 and
# 01132022), or the two pooled runs combined (id should be "pooled")
sample_id = sys.argv[1]

input_file = "../celltypist/%s_input_data.csv" % sample_id

# Run celltypist
predictions = celltypist.annotate(input_file, majority_voting = True)

# Save output files
label_file = "../celltypist/%s_predicted_labels.csv" % sample_id
predictions.predicted_labels.to_csv(label_file, sep="\t")
prob_file = "../celltypist/%s_probability_matrix.csv" % sample_id
predictions.probability_matrix.to_csv(prob_file, sep="\t")
