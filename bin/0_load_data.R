#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
path = args[1]

###############################
###   0 - Getting Started   ###
###############################

##############################
###   0.1 - Loading Data   ###
##############################

library(GeomxTools)

datadir <- system.file("extdata", "WTA_NGS_Example", package="GeoMxWorkflows")

DCCFiles <- dir(
  file.path(path, "input_data", "dccs"),
  pattern = ".dcc$",
  full.names = TRUE,
  recursive = TRUE
)
PKCFiles <- file.path(path, "input_data", "Hsa_WTA_v1.0.pkc")
SampleAnnotationFile <- file.path(path, "input_data", "kidney_AOI_Annotations_all_vignette.xlsx")

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
save.image(sprintf('%s/image/0_load_data.RData', path))
