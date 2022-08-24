#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
pathBase = args[1]
minLOQ = as.integer(args[2])
cutoffLOQ = as.integer(args[3])

load(sprintf('%s/image/2_3_aggregate_counts.RData', pathBase))

#################################################
###   Section 2.4 - Limit of Quantification   ###
#################################################

library(Biobase)


# Calculate LOQ per module tested
LOQ <- data.frame(row.names = colnames(target_data))
for(module in modules) {
  vars <- paste0(c("NegGeoMean_", "NegGeoSD_"),
                 module)
  if(all(vars[1:2] %in% colnames(pData(target_data)))) {
    LOQ[, module] <-
      pmax(minLOQ,
           pData(target_data)[, vars[1]] * 
             pData(target_data)[, vars[2]] ^ cutoffLOQ)
  }
}
pData(target_data)$LOQ <- LOQ

# Save image
save.image(sprintf('%s/image/2_4_limit_of_quantification.RData', pathBase))
