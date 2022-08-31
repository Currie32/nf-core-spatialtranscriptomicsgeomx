#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
pathBase = args[1]
pathInputData = args[2]
pkcFiles = args[3]
annotationFile = args[4]

###############################
###   1 - Getting Started   ###
###############################

##############################
###   1.1 - Loading Data   ###
##############################

library(GeomxTools)

DCCFiles <- dir(
  file.path(pathBase, pathInputData, "dccs"),
  pattern = ".dcc$",
  full.names = TRUE,
  recursive = TRUE
)
PKCFiles <- file.path(pathBase, pathInputData, pkcFiles)
SampleAnnotationFile <- file.path(pathBase, pathInputData, annotationFile)

data <- readNanoStringGeoMxSet(
  dccFiles = DCCFiles,
  pkcFiles = PKCFiles,
  phenoDataFile = SampleAnnotationFile,
  phenoDataSheet = "Template",
  phenoDataDccColName = "Sample_ID",
  protocolDataColNames = c("aoi", "roi"),
  experimentDataColNames = c("panel")
)

# Save image
save.image(sprintf('%s/image/1_load_data.RData', pathBase))
