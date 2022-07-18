#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
path = args[1]

###############################
###   2 - Getting Started   ###
###############################

##############################
###   2.1 - Loading Data   ###
##############################

library(GeomxTools)

datadir <- system.file("extdata", "WTA_NGS_Example", package="GeoMxWorkflows")

DCCFiles <- dir(
  file.path(datadir, "dccs"),
  pattern = ".dcc$",
  full.names = TRUE,
  recursive = TRUE
)

PKCFiles <- unzip(zipfile = dir(
  file.path(datadir, "pkcs"),
  pattern = ".zip$",
  full.names = TRUE,
  recursive = TRUE
))

SampleAnnotationFile <- dir(
  file.path(datadir, "annotation"),
  pattern = ".xlsx$",
  full.names = TRUE,
  recursive = TRUE
)

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
save.image(sprintf('%s/image/2_load_data.RData', path))
