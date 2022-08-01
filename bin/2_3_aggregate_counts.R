#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
path = args[1]

load(sprintf('%s/image/2_2_probe_qc.RData', path))

######################################################
###   Section 2.3 - Create Gene-level Count Data   ###
######################################################

library(GeomxTools)

# collapse to targets
target_data <- aggregateCounts(data)

# Save the dimensions of target_data
write.csv(
  data.frame(dim(target_data)),
  sprintf("%s/data/2_3_dimensions_target_data.csv", path),
  row.names=FALSE
)

# Save image
save.image(sprintf('%s/image/2_3_aggregate_counts.RData', path))
