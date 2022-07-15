#!/usr/bin/env Rscript

load('~/Imperial/nf-core-spatialtranscriptomicsgeomx/data/3_sample_overview.RData')

############################
###   4.1 - Segment QC   ###
############################

#####################################
###   4.1.1 - Select Segment QC   ###
#####################################

library(GeomxTools)

# Shift counts to one
demoData <- shiftCountsOne(demoData, useDALogic = TRUE)

# Default QC cutoffs are commented in () adjacent to the respective parameters
# study-specific values were selected after visualizing the QC results in more
# detail below
QC_params <-
  list(minSegmentReads = 1000, # Minimum number of reads (1000)
       percentTrimmed = 80,    # Minimum % of reads trimmed (80%)
       percentStitched = 80,   # Minimum % of reads stitched (80%)
       percentAligned = 75,    # Minimum % of reads aligned (80%)
       percentSaturation = 50, # Minimum sequencing saturation (50%)
       minNegativeCount = 1,   # Minimum negative control counts (10)
       maxNTCCount = 9000,     # Maximum counts observed in NTC well (1000)
       minNuclei = 20,         # Minimum # of nuclei estimated (100)
       minArea = 1000)         # Minimum segment area (5000)
demoData <-
  setSegmentQCFlags(demoData, 
                    qcCutoffs = QC_params)        

# Collate QC Results
QCResults <- protocolData(demoData)[["QCFlags"]]
flag_columns <- colnames(QCResults)
QC_Summary <- data.frame(Pass = colSums(!QCResults[, flag_columns]),
                         Warning = colSums(QCResults[, flag_columns]))
QCResults$QCStatus <- apply(QCResults, 1L, function(x) {
  ifelse(sum(x) == 0L, "PASS", "WARNING")
})
QC_Summary["TOTAL FLAGS", ] <-
  c(sum(QCResults[, "QCStatus"] == "PASS"),
    sum(QCResults[, "QCStatus"] == "WARNING"))


################################################
###   Section 4.1.2 - Visualize Segment QC   ###
################################################

library(ggplot2)
library(knitr)

# Graphical summaries of QC statistics plot function
QC_histogram <- function(assay_data = NULL,
                         annotation = NULL,
                         fill_by = NULL,
                         thr = NULL,
                         scale_trans = NULL) {
  plt <- ggplot(assay_data,
                aes_string(x = paste0("unlist(`", annotation, "`)"),
                           fill = fill_by)) +
    geom_histogram(bins = 50) +
    geom_vline(xintercept = thr, lty = "dashed", color = "black") +
    theme_bw() + guides(fill = "none") +
    facet_wrap(as.formula(paste("~", fill_by)), nrow = 4) +
    labs(x = annotation, y = "Segments, #", title = annotation)
  if(!is.null(scale_trans)) {
    plt <- plt +
      scale_x_continuous(trans = scale_trans)
  }
  plt
}

# Column to plot
col_by <- "segment"

QC_histogram(sData(demoData), "Trimmed (%)", col_by, 80)
ggsave("~/Imperial/nf-core-spatialtranscriptomicsgeomx/plots/4_1_2_trimmed_percentage.png", device='png')

QC_histogram(sData(demoData), "Stitched (%)", col_by, 80)
ggsave("~/Imperial/nf-core-spatialtranscriptomicsgeomx/plots/4_1_2_stitched_percentage.png", device='png')

QC_histogram(sData(demoData), "Aligned (%)", col_by, 75)
ggsave("~/Imperial/nf-core-spatialtranscriptomicsgeomx/plots/4_1_2_aligned_percentage.png", device='png')

QC_histogram(sData(demoData), "Saturated (%)", col_by, 50) +
  labs(title = "Sequencing Saturation (%)",
       x = "Sequencing Saturation (%)")
ggsave("~/Imperial/nf-core-spatialtranscriptomicsgeomx/plots/4_1_2_satured_percentage.png", device='png')

QC_histogram(sData(demoData), "area", col_by, 1000, scale_trans = "log10")
ggsave("~/Imperial/nf-core-spatialtranscriptomicsgeomx/plots/4_1_2_area.png", device='png')

QC_histogram(sData(demoData), "nuclei", col_by, 20)
ggsave("~/Imperial/nf-core-spatialtranscriptomicsgeomx/plots/4_1_2_nuclei.png", device='png')


# calculate the negative geometric means for each module
negativeGeoMeans <- 
  esBy(negativeControlSubset(demoData), 
       GROUP = "Module", 
       FUN = function(x) { 
         assayDataApply(x, MARGIN = 2, FUN = ngeoMean, elt = "exprs") 
       }) 
protocolData(demoData)[["NegGeoMean"]] <- negativeGeoMeans

# explicitly copy the Negative geoMeans from sData to pData
negCols <- paste0("NegGeoMean_", modules)
pData(demoData)[, negCols] <- sData(demoData)[["NegGeoMean"]]
for(ann in negCols) {
  plt <- QC_histogram(pData(demoData), ann, col_by, 2, scale_trans = "log10")
  print(plt)
  ggsave(sprintf("~/Imperial/nf-core-spatialtranscriptomicsgeomx/plots/4_1_2_NegGeoMean_%s.png", ann), device='png')
}

# detatch neg_geomean columns ahead of aggregateCounts call
pData(demoData) <- pData(demoData)[, !colnames(pData(demoData)) %in% negCols]

# show all NTC values, Freq = # of Segments with a given NTC count:
ntc_count_table <- kable(
  table(NTC_Count = sData(demoData)$NTC),
  col.names = c("NTC Count", "# of Segments"),
  caption="Number of Segments with a given NTC count"
)
file_conn <- file("~/Imperial/nf-core-spatialtranscriptomicsgeomx/data/4_1_2_table_ntc_count.txt")
writeLines(ntc_count_table, file_conn)
close(file_conn)

qc_summary_table <- kable(QC_Summary, caption="QC Summary Table for each Segment")
file_conn <- file("~/Imperial/nf-core-spatialtranscriptomicsgeomx/data/4_1_2_table_qc_summary.txt")
writeLines(qc_summary_table, file_conn)
close(file_conn)

###################################################
###   Section 4.1.3 - Remove flagged segments   ###
###################################################

demoData <- demoData[, QCResults$QCStatus == "PASS"]
write.csv(
  data.frame(dim(demoData)),
  "~/Imperial/nf-core-spatialtranscriptomicsgeomx/data/4_1_3_dimensions_after_segment_qc.csv",
)

# Save image
save.image('~/Imperial/nf-core-spatialtranscriptomicsgeomx/image/4_1_segment_qc.RData')

