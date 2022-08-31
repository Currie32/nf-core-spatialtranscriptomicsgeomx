#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
pathBase = args[1]

load(sprintf('%s/image/3_2_probe_qc.RData', pathBase))

######################################################
###   Section 3.3 - Create Gene-level Count Data   ###
######################################################

library(GeomxTools)

# collapse to targets
target_data <- aggregateCounts(data)

# Save the dimensions of target_data
write.csv(
  data.frame(dim(target_data)),
  sprintf("%s/data/3_3_dimensions_target_data.csv", pathBase),
  row.names=FALSE
)

# Save image
save.image(sprintf('%s/image/3_3_aggregate_counts.RData', pathBase))
