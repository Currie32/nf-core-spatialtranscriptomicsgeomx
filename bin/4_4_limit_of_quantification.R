#!/usr/bin/env Rscript

load('~/Imperial/nf-core-spatialtranscriptomicsgeomx/image/4_3_aggregate_counts.RData')

#################################################
###   Section 4.4 - Limit of Quantification   ###
#################################################

library(Biobase)

# Define LOQ SD threshold and minimum value
cutoff <- 2
minLOQ <- 2

# Calculate LOQ per module tested
LOQ <- data.frame(row.names = colnames(target_demoData))
for(module in modules) {
  vars <- paste0(c("NegGeoMean_", "NegGeoSD_"),
                 module)
  if(all(vars[1:2] %in% colnames(pData(target_demoData)))) {
    LOQ[, module] <-
      pmax(minLOQ,
           pData(target_demoData)[, vars[1]] * 
             pData(target_demoData)[, vars[2]] ^ cutoff)
  }
}
pData(target_demoData)$LOQ <- LOQ

# Save image
save.image('~/Imperial/nf-core-spatialtranscriptomicsgeomx/image/4_4_limit_of_quantification.RData')
