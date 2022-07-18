#!/usr/bin/env Rscript

load('~/Imperial/nf-core-spatialtranscriptomicsgeomx/image/4_2_probe_qc.RData')

######################################################
###   Section 4.3 - Create Gene-level Count Data   ###
######################################################

library(GeomxTools)

# collapse to targets
target_data <- aggregateCounts(data)

# Save the dimensions of target_data
write.csv(
  data.frame(dim(target_data)),
  "~/Imperial/nf-core-spatialtranscriptomicsgeomx/data/dimensions_target_data.csv",
)

# Save image
save.image('~/Imperial/nf-core-spatialtranscriptomicsgeomx/image/4_3_aggregate_counts.RData')
