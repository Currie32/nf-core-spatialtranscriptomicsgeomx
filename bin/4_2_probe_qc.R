#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
path = args[1]

load(sprintf('%s/image/4_1_segment_qc.RData', path))

##################################
###   Section 4.2 - Probe QC   ###
##################################

##############################################
###   Section 4.2.1 - Set Probe QC Flags   ###
##############################################

library(Biobase)
library(GeomxTools)

# Generally keep the qcCutoffs parameters unchanged. Set removeLocalOutliers to 
# FALSE if you do not want to remove local outliers
data <- setBioProbeQCFlags(data,
                           qcCutoffs = list(minProbeRatio = 0.1,
                                            percentFailGrubbs = 20),
                           removeLocalOutliers = TRUE)

ProbeQCResults <- fData(data)[["QCFlags"]]

# Define QC table for Probe QC
qc_df <- data.frame(Passed = sum(rowSums(ProbeQCResults[, -1]) == 0),
                    Global = sum(ProbeQCResults$GlobalGrubbsOutlier),
                    Local = sum(rowSums(ProbeQCResults[, -2:-1]) > 0
                                & !ProbeQCResults$GlobalGrubbsOutlier))


##################################################
###   Section 4.2.2 - Exclude Outlier Probes   ###
##################################################

#Subset object to exclude all that did not pass Ratio & Global testing
ProbeQCPassed <- 
  subset(data, 
         fData(data)[["QCFlags"]][,c("LowProbeRatio")] == FALSE &
           fData(data)[["QCFlags"]][,c("GlobalGrubbsOutlier")] == FALSE)

# Save the dimensions of ProbeQCPassed
write.csv(
  data.frame(dim(ProbeQCPassed)),
  sprintf("%s/data/4_2_2_dimensions_after_probe_qc.csv", path)
)

data <- ProbeQCPassed 

# Save image
save.image(sprintf('%s/image/4_2_probe_qc.RData', path))
