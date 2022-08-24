#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
pathBase = args[1]
minProbeRatio = as.numeric(args[2])
percentFailGrubbs = as.integer(args[3])

load(sprintf('%s/image/2_1_segment_qc.RData', pathBase))

##################################
###   Section 2.2 - Probe QC   ###
##################################

##############################################
###   Section 2.2.1 - Set Probe QC Flags   ###
##############################################

library(Biobase)
library(GeomxTools)

# Generally keep the qcCutoffs parameters unchanged. Set removeLocalOutliers to 
# FALSE if you do not want to remove local outliers
data <- setBioProbeQCFlags(data,
                           qcCutoffs = list(minProbeRatio = minProbeRatio,
                                            percentFailGrubbs = percentFailGrubbs),
                           removeLocalOutliers = TRUE)

ProbeQCResults <- fData(data)[["QCFlags"]]

# Define QC table for Probe QC
qc_df <- data.frame(Passed = sum(rowSums(ProbeQCResults[, -1]) == 0),
                    Global = sum(ProbeQCResults$GlobalGrubbsOutlier),
                    Local = sum(rowSums(ProbeQCResults[, -2:-1]) > 0
                                & !ProbeQCResults$GlobalGrubbsOutlier))


##################################################
###   Section 2.2.2 - Exclude Outlier Probes   ###
##################################################

#Subset object to exclude all that did not pass Ratio & Global testing
ProbeQCPassed <- 
  subset(data, 
         fData(data)[["QCFlags"]][,c("LowProbeRatio")] == FALSE &
           fData(data)[["QCFlags"]][,c("GlobalGrubbsOutlier")] == FALSE)

# Save the dimensions of ProbeQCPassed
write.csv(
  data.frame(dim(ProbeQCPassed)),
  sprintf("%s/data/2_2_2_dimensions_after_probe_qc.csv", pathBase)
)

data <- ProbeQCPassed 

# Save image
save.image(sprintf('%s/image/2_2_probe_qc.RData', pathBase))
