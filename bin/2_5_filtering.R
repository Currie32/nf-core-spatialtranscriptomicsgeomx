#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
path = args[1]

load(sprintf('%s/image/2_4_limit_of_quantification.RData', path))

###################################
###   Section 2.5 - Filtering   ###
###################################

library(Biobase)

LOQ_Mat <- c()
for(module in modules) {
  ind <- fData(target_data)$Module == module
  Mat_i <- t(esApply(target_data[ind, ], MARGIN = 1,
                     FUN = function(x) {
                       x > LOQ[, module]
                     }))
  LOQ_Mat <- rbind(LOQ_Mat, Mat_i)
}
# ensure ordering since this is stored outside of the geomxSet
LOQ_Mat <- LOQ_Mat[fData(target_data)$TargetName, ]


##################################################
###   Section 2.5.1 - Segment Gene Detection   ###
##################################################

library(dplyr)
library(ggforce)
library(knitr)
library(tidyr)

# Save detection rate information to pheno data
pData(target_data)$GenesDetected <- 
  colSums(LOQ_Mat, na.rm = TRUE)
pData(target_data)$GeneDetectionRate <-
  pData(target_data)$GenesDetected / nrow(target_data)

# Determine detection thresholds: 1%, 5%, 10%, 15%, >15%
pData(target_data)$DetectionThreshold <- 
  cut(pData(target_data)$GeneDetectionRate,
      breaks = c(0, 0.01, 0.05, 0.1, 0.15, 1),
      labels = c("<1%", "1-5%", "5-10%", "10-15%", ">15%"))

# stacked bar plot of different cut points (1%, 5%, 10%, 15%)
ggplot(pData(target_data),
       aes(x = DetectionThreshold)) +
  geom_bar(aes(fill = region)) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Gene Detection Rate",
       y = "Segments, #",
       fill = "Segment Type")

ggsave(sprintf("%s/plots/2_5_1_gene_detection_rate_by_segment.png", path), device='png')


# cut percent genes detected at 1, 5, 10, 15
gene_detection_rates <- data.frame(
  table(pData(target_data)$DetectionThreshold,
        pData(target_data)$class)
)
gene_detection_rates <- gene_detection_rates %>%
  pivot_wider(id_cols="Var1", names_from = "Var2", values_from = "Freq")

names(gene_detection_rates)[names(gene_detection_rates) == 'Var1'] <- 'Gene detection rate'

write.csv(
  gene_detection_rates,
  sprintf("%s/data/2_5_1_table_gene_detection_rate_by_kidney_tissue_type.csv", path),
  row.names=FALSE
)


target_data <- target_data[, pData(target_data)$GeneDetectionRate >= .1]

# Save the dimensions of target_data
write.csv(
  data.frame(dim(target_data)),
  sprintf("%s/data/2_5_1_dimensions_target_data_after_gene_detection_rate_filter.csv", path)
)

# select the annotations we want to show, use `` to surround column names with
# spaces or special symbols
count_mat <- count(pData(data), `slide name`, class, region, segment)
# simplify the slide names
count_mat$`slide name` <- 
  gsub("disease", "d",
       gsub("normal", "n", count_mat$`slide name`))
# gather the data and plot in order: class, slide name, region, segment
test_gr <- gather_set_data(count_mat, 1:4)
test_gr$x <-
  factor(test_gr$x,
         levels = c("class", "slide name", "region", "segment"))
# plot Sankey
ggplot(test_gr, aes(x, id = id, split = y, value = n)) +
  geom_parallel_sets(aes(fill = region), alpha = 0.5, axis.width = 0.1) +
  geom_parallel_sets_axes(axis.width = 0.2) +
  geom_parallel_sets_labels(color = "white", size = 5) +
  theme_classic(base_size = 17) + 
  theme(legend.position = "bottom",
        axis.ticks.y = element_blank(),
        axis.line = element_blank(),
        axis.text.y = element_blank()) +
  scale_y_continuous(expand = expansion(0)) + 
  scale_x_discrete(expand = expansion(0)) +
  labs(x = "", y = "") +
  annotate(geom = "segment", x = 4.25, xend = 4.25, y = 20, 
           yend = 120, lwd = 2) +
  annotate(geom = "text", x = 4.19, y = 70, angle = 90, size = 5,
           hjust = 0.5, label = "100 segments")

ggsave(sprintf("%s/plots/2_5_1_sample_overview_sankey_after_gene_detection_rate_filter.png", path), device='png')


###############################################
###   Section 2.5.2 - Gene Detection Rate   ###
###############################################

library(scales) # for percent

# Calculate detection rate:
LOQ_Mat <- LOQ_Mat[, colnames(target_data)]
fData(target_data)$DetectedSegments <- rowSums(LOQ_Mat, na.rm = TRUE)
fData(target_data)$DetectionRate <-
  fData(target_data)$DetectedSegments / nrow(pData(target_data))

# Gene of interest detection table
goi <- c("PDCD1", "CD274", "IFNG", "CD8A", "CD68", "EPCAM",
         "KRT18", "NPHS1", "NPHS2", "CALB1", "CLDN8")
goi_df <- data.frame(
  Gene = goi,
  Number = fData(target_data)[goi, "DetectedSegments"],
  DetectionRate = percent(fData(target_data)[goi, "DetectionRate"]))

write.csv(
  goi_df,
  sprintf("%s/data/2_5_2_detection_rate_for_genes_of_interest.csv", path),
  row.names=FALSE
)


##########################################
###   Section 2.5.3 - Gene Filtering   ###
##########################################

# Plot detection rate:
plot_detect <- data.frame(Freq = c(1, 5, 10, 20, 30, 50))
plot_detect$Number <-
  unlist(lapply(c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5),
                function(x) {sum(fData(target_data)$DetectionRate >= x)}))
plot_detect$Rate <- plot_detect$Number / nrow(fData(target_data))
rownames(plot_detect) <- plot_detect$Freq

ggplot(plot_detect, aes(x = as.factor(Freq), y = Rate, fill = Rate)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = formatC(Number, format = "d", big.mark = ",")),
            vjust = 1.6, color = "black", size = 4) +
  scale_fill_gradient2(low = "orange2", mid = "lightblue",
                       high = "dodgerblue3", midpoint = 0.65,
                       limits = c(0,1),
                       labels = scales::percent) +
  theme_bw() +
  scale_y_continuous(labels = scales::percent, limits = c(0,1),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = "% of Segments",
       y = "Genes Detected, % of Panel > LOQ")

ggsave(sprintf("%s/plots/2_5_3_genes_detected_per_percentage_of_segments.png", path), device='png')


# Subset to target genes detected in at least 10% of the samples.
#   Also manually include the negative control probe, for downstream use
negativeProbefData <- subset(fData(target_data), CodeClass == "Negative")
neg_probes <- unique(negativeProbefData$TargetName)
target_data <- 
  target_data[fData(target_data)$DetectionRate >= 0.1 |
                fData(target_data)$TargetName %in% neg_probes, ]

# Save the dimensions of target_data
write.csv(
  data.frame(dim(target_data)),
  sprintf("%s/data/2_5_3_dimensions_target_data_after_gene_detection_filter.csv", path)
)

# retain only detected genes of interest
goi <- goi[goi %in% rownames(target_data)]

# Save image
save.image(sprintf('%s/image/2_5_filtering.RData', path))
