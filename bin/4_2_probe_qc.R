#!/usr/bin/env Rscript

load('~/Imperial/nf-core-spatialtranscriptomicsgeomx/image/4_1_segment_qc.RData')

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
demoData <- setBioProbeQCFlags(demoData, 
                               qcCutoffs = list(minProbeRatio = 0.1,
                                                percentFailGrubbs = 20), 
                               removeLocalOutliers = TRUE)

ProbeQCResults <- fData(demoData)[["QCFlags"]]

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
  subset(demoData, 
         fData(demoData)[["QCFlags"]][,c("LowProbeRatio")] == FALSE &
           fData(demoData)[["QCFlags"]][,c("GlobalGrubbsOutlier")] == FALSE)

# Save the dimensions of ProbeQCPassed
write.csv(
  data.frame(dim(ProbeQCPassed)),
  "~/Imperial/nf-core-spatialtranscriptomicsgeomx/data/4_2_2_dimensions_after_probe_qc.csv",
)

demoData <- ProbeQCPassed 

# Save image
save.image('~/Imperial/nf-core-spatialtranscriptomicsgeomx/image/4_2_probe_qc.RData')
