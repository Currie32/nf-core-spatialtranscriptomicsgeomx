#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
path = args[1]

load(sprintf('%s/image/3_sample_overview.RData', path))

############################
###   4.1 - Segment QC   ###
############################

#####################################
###   4.1.1 - Select Segment QC   ###
#####################################

library(GeomxTools)

# Shift counts to one
data <- shiftCountsOne(data, useDALogic = TRUE)

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
data <-
  setSegmentQCFlags(data, 
                    qcCutoffs = QC_params)        

# Collate QC Results
QCResults <- protocolData(data)[["QCFlags"]]
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

QC_histogram(sData(data), "Trimmed (%)", col_by, 80)
ggsave(sprintf("%s/plots/4_1_2_trimmed_percentage.png", path), device='png')

QC_histogram(sData(data), "Stitched (%)", col_by, 80)
ggsave(sprintf("%s/plots/4_1_2_stitched_percentage.png", path), device='png')

QC_histogram(sData(data), "Aligned (%)", col_by, 75)
ggsave(sprintf("%s/plots/4_1_2_aligned_percentage.png", path), device='png')

QC_histogram(sData(data), "Saturated (%)", col_by, 50) +
  labs(title = "Sequencing Saturation (%)",
       x = "Sequencing Saturation (%)")
ggsave(sprintf("%s/plots/4_1_2_satured_percentage.png", path), device='png')

QC_histogram(sData(data), "area", col_by, 1000, scale_trans = "log10")
ggsave(sprintf("%s/plots/4_1_2_area.png", path), device='png')

QC_histogram(sData(data), "nuclei", col_by, 20)
ggsave(sprintf("%s/plots/4_1_2_nuclei.png", path), device='png')


# calculate the negative geometric means for each module
negativeGeoMeans <- 
  esBy(negativeControlSubset(data), 
       GROUP = "Module", 
       FUN = function(x) { 
         assayDataApply(x, MARGIN = 2, FUN = ngeoMean, elt = "exprs") 
       }) 
protocolData(data)[["NegGeoMean"]] <- negativeGeoMeans

# explicitly copy the Negative geoMeans from sData to pData
negCols <- paste0("NegGeoMean_", modules)
pData(data)[, negCols] <- sData(data)[["NegGeoMean"]]
for(ann in negCols) {
  plt <- QC_histogram(pData(data), ann, col_by, 2, scale_trans = "log10")
  print(plt)
  ggsave(sprintf("%s/plots/4_1_2_NegGeoMean_%s.png", path, ann), device='png')
}

# detatch neg_geomean columns ahead of aggregateCounts call
pData(data) <- pData(data)[, !colnames(pData(data)) %in% negCols]

# show all NTC values, Freq = # of Segments with a given NTC count:
ntc_count_table <- kable(
  table(NTC_Count = sData(data)$NTC),
  col.names = c("NTC Count", "# of Segments"),
  caption="Number of Segments with a given NTC count"
)
file_conn <- file(sprintf("%s/data/4_1_2_table_ntc_count.txt", path))
writeLines(ntc_count_table, file_conn)
close(file_conn)

qc_summary_table <- kable(QC_Summary, caption="QC Summary Table for each Segment")
file_conn <- file(sprintf("%s/data/4_1_2_table_qc_summary.txt", path))
writeLines(qc_summary_table, file_conn)
close(file_conn)

###################################################
###   Section 4.1.3 - Remove flagged segments   ###
###################################################

data <- data[, QCResults$QCStatus == "PASS"]
write.csv(
  data.frame(dim(data)),
  sprintf("%s/data/4_1_3_dimensions_after_segment_qc.csv", path)
)

# Save image
save.image(sprintf('%s/image/4_1_segment_qc.RData', path))
