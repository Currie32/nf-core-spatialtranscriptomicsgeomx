#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
path = args[1]

load(sprintf('%s/image/7_differential_expression.RData', path))

############################################
###   Section 8 - Visualizing DE Genes   ###
############################################

#######################################
###   Section 8.1 - Volcano Plots   ###
#######################################

library(ggrepel)

# Categorize Results based on P-value & FDR for plotting
results$Color <- "NS or FC < 0.5"
results$Color[results$`Pr(>|t|)` < 0.05] <- "P < 0.05"
results$Color[results$FDR < 0.05] <- "FDR < 0.05"
results$Color[results$FDR < 0.001] <- "FDR < 0.001"
results$Color[abs(results$Estimate) < 0.5] <- "NS or FC < 0.5"
results$Color <- factor(results$Color,
                        levels = c("NS or FC < 0.5", "P < 0.05",
                                   "FDR < 0.05", "FDR < 0.001"))

# pick top genes for either side of volcano to label
# order genes for convenience:
results$invert_P <- (-log10(results$`Pr(>|t|)`)) * sign(results$Estimate)

top_g <- c()
for(cond in c("DKD", "normal")) {
  ind <- results$Subset == cond
  top_g <- c(top_g,
             results[ind, 'Gene'][
               order(results[ind, 'invert_P'], decreasing = TRUE)[1:15]],
             results[ind, 'Gene'][
               order(results[ind, 'invert_P'], decreasing = FALSE)[1:15]])
}

top_g <- unique(top_g)
results <- results[, -1*ncol(results)] # remove invert_P from matrix

# Graph results
ggplot(results,
       aes(x = Estimate, y = -log10(`Pr(>|t|)`),
           color = Color, label = Gene)) +
  geom_vline(xintercept = c(0.5, -0.5), lty = "dashed") +
  geom_hline(yintercept = -log10(0.05), lty = "dashed") +
  geom_point() +
  labs(x = "Enriched in Tubules <- log2(FC) -> Enriched in Glomeruli",
       y = "Significance, -log10(P)",
       color = "Significance") +
  scale_color_manual(values = c(`FDR < 0.001` = "dodgerblue",
                                `FDR < 0.05` = "lightblue",
                                `P < 0.05` = "orange2",
                                `NS or FC < 0.5` = "gray"),
                     guide = guide_legend(override.aes = list(size = 4))) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  geom_text_repel(data = subset(results, Gene %in% top_g & FDR < 0.001),
                  size = 4, point.padding = 0.15, color = "black",
                  min.segment.length = .1, box.padding = .2, lwd = 2,
                  max.overlaps = 50) +
  theme_bw(base_size = 16) +
  theme(legend.position = "bottom") +
  facet_wrap(~Subset, scales = "free_y")

ggsave(sprintf("%s/plots/8_1_volcano_plots.png", path), device='png')


####################################################
###   Section 8.2 - Plotting Genes of Interest   ###
####################################################

library(Biobase)
library(knitr)

table <- kable(subset(results, Gene %in% c('PDHA1','ITGB1')), row.names = FALSE)

file_conn <- file(sprintf("%s/data/8_2_model_results_genes_of_interest.txt", path))
writeLines(table, file_conn)
close(file_conn)

# show expression for a single target: PDHA1
ggplot(pData(target_data),
       aes(x = region, fill = region,
           y = assayDataElement(target_data["PDHA1", ],
                                elt = "q_norm"))) +
  geom_violin() +
  geom_jitter(width = .2) +
  labs(y = "PDHA1 Expression") +
  scale_y_continuous(trans = "log2") +
  facet_wrap(~class) +
  theme_bw()

ggsave(sprintf("%s/plots/8_2_violin_plot_gene_expression.png", path), device='png')


glom <- pData(target_data)$region == "glomerulus"

# show expression of PDHA1 vs ITGB1
ggplot(pData(target_data),
       aes(x = assayDataElement(target_data["PDHA1", ],
                                elt = "q_norm"),
           y = assayDataElement(target_data["ITGB1", ],
                                elt = "q_norm"),
           color = region)) +
  geom_vline(xintercept =
               max(assayDataElement(target_data["PDHA1", glom],
                                    elt = "q_norm")),
             lty = "dashed", col = "darkgray") +
  geom_hline(yintercept =
               max(assayDataElement(target_data["ITGB1", !glom],
                                    elt = "q_norm")),
             lty = "dashed", col = "darkgray") +
  geom_point(size = 3) +
  theme_bw() +
  scale_x_continuous(trans = "log2") + 
  scale_y_continuous(trans = "log2") +
  labs(x = "PDHA1 Expression", y = "ITGB1 Expression") +
  facet_wrap(~class)

ggsave(sprintf("%s/plots/8_2_expression_patterns.png", path), device='png')


######################################################
###   Section 8.3 - Heatmap of Significant Genes   ###
######################################################

library(pheatmap)

# select top significant genes based on significance, plot with pheatmap
GOI <- unique(subset(results, `FDR` < 0.001)$Gene)
pheatmap(log2(assayDataElement(target_data[GOI, ], elt = "q_norm")),
         scale = "row", 
         show_rownames = FALSE, show_colnames = FALSE,
         border_color = NA,
         clustering_method = "average",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         cutree_cols = 2, cutree_rows = 2,
         breaks = seq(-3, 3, 0.05),
         color = colorRampPalette(c("purple3", "black", "yellow2"))(120),
         annotation_col = pData(target_data)[, c("region", "class")])

ggsave(sprintf("%s/plots/8_3_heatmap_significant_genes.png", path), device='png')


###################################
###   Section 9.3.1 - MA Plot   ###
###################################

results$MeanExp <- rowMeans(assayDataElement(target_data, elt = "q_norm"))

top_g2 <- results$Gene[results$Gene %in% top_g &
                         results$FDR < 0.001 &
                         abs(results$Estimate) > .5 &
                         results$MeanExp > quantile(results$MeanExp, 0.9)]

ggplot(subset(results, !Gene %in% neg_probes),
       aes(x = MeanExp, y = Estimate,
           size = -log10(`Pr(>|t|)`),
           color = Color, label = Gene)) +
  geom_hline(yintercept = c(0.5, -0.5), lty = "dashed") +
  scale_x_continuous(trans = "log2") +
  geom_point(alpha = 0.5) + 
  labs(y = "Enriched in Glomeruli <- log2(FC) -> Enriched in Tubules",
       x = "Mean Expression",
       color = "Significance") +
  scale_color_manual(values = c(`FDR < 0.001` = "dodgerblue",
                                `FDR < 0.05` = "lightblue",
                                `P < 0.05` = "orange2",
                                `NS or FC < 0.5` = "gray")) +
  geom_text_repel(data = subset(results, Gene %in% top_g2),
                  size = 4, point.padding = 0.15, color = "black",
                  min.segment.length = .1, box.padding = .2, lwd = 2) +
  theme_bw(base_size = 16) +
  facet_wrap(~Subset, nrow = 2, ncol = 1)

# Saved under section 8.4 instead of 9.3.1 to avoid confusion of which script produces this plot
ggsave(sprintf("%s/plots/8_4_ma_plot.png", path), device='png')

# Save image
save.image(sprintf('%s/image/8_visualising_de_genes.RData', path))
