#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
path = args[1]

load(sprintf('%s/image/4_5_filtering.RData', path))

#####################################
###   Section 5 - Normalization   ###
#####################################

library(Biobase)
library(cowplot)   # for plot_grid
library(reshape2)  # for melt

# Graph Q3 value vs negGeoMean of Negatives
ann_of_interest <- "region"
Stat_data <- 
  data.frame(row.names = colnames(exprs(target_data)),
             Segment = colnames(exprs(target_data)),
             Annotation = pData(target_data)[, ann_of_interest],
             Q3 = unlist(apply(exprs(target_data), 2,
                               quantile, 0.75, na.rm = TRUE)),
             NegProbe = exprs(target_data)[neg_probes, ])
Stat_data_m <- melt(Stat_data, measure.vars = c("Q3", "NegProbe"),
                    variable.name = "Statistic", value.name = "Value")

plt1 <- ggplot(Stat_data_m,
               aes(x = Value, fill = Statistic)) +
  geom_histogram(bins = 40) + theme_bw() +
  scale_x_continuous(trans = "log2") +
  facet_wrap(~Annotation, nrow = 1) + 
  scale_fill_brewer(palette = 3, type = "qual") +
  labs(x = "Counts", y = "Segments, #")

plt2 <- ggplot(Stat_data,
               aes(x = NegProbe, y = Q3, color = Annotation)) +
  geom_abline(intercept = 0, slope = 1, lty = "dashed", color = "darkgray") +
  geom_point() + guides(color = "none") + theme_bw() +
  scale_x_continuous(trans = "log2") + 
  scale_y_continuous(trans = "log2") +
  theme(aspect.ratio = 1) +
  labs(x = "Negative Probe GeoMean, Counts", y = "Q3 Value, Counts")

plt3 <- ggplot(Stat_data,
               aes(x = NegProbe, y = Q3 / NegProbe, color = Annotation)) +
  geom_hline(yintercept = 1, lty = "dashed", color = "darkgray") +
  geom_point() + theme_bw() +
  scale_x_continuous(trans = "log2") + 
  scale_y_continuous(trans = "log2") +
  theme(aspect.ratio = 1) +
  labs(x = "Negative Probe GeoMean, Counts", y = "Q3/NegProbe Value, Counts")

btm_row <- plot_grid(plt2, plt3, nrow = 1, labels = c("B", ""),
                     rel_widths = c(0.43,0.57))
plot_grid(plt1, btm_row, ncol = 1, labels = c("A", ""))

ggsave(sprintf("%s/plots/5_q3_values_vs_neg_geo_mean_of_neg_probes.png", path), device='png')


# Q3 norm (75th percentile) for WTA/CTA  with or without custom spike-ins
target_data <- normalize(target_data,
                         norm_method = "quant", 
                         desiredQuantile = .75,
                         toElt = "q_norm")

# Background normalization for WTA/CTA without custom spike-in
target_data <- normalize(target_data,
                         norm_method = "neg", 
                         fromElt = "exprs",
                         toElt = "neg_norm")

# visualize the first 10 segments with each normalization method
boxplot(exprs(target_data)[,1:10],
        col = "#9EDAE5", main = "Raw Counts",
        log = "y", names = 1:10, xlab = "Segment",
        ylab = "Counts, Raw")

ggsave(sprintf("%s/plots/5_segment_counts_raw.png", path), device='png')


boxplot(assayDataElement(target_data[,1:10], elt = "q_norm"),
        col = "#2CA02C", main = "Q3 Norm Counts",
        log = "y", names = 1:10, xlab = "Segment",
        ylab = "Counts, Q3 Normalized")

ggsave(sprintf("%s/plots/5_segment_counts_normalised.png", path), device='png')

# Save image
save.image(sprintf('%s/image/5_normalisation.RData', path))
