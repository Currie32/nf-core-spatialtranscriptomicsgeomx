#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
path = args[1]

load(sprintf('%s/image/3_normalisation.RData', path))

#############################################
###   Section 4 - Unsupervised Analysis   ###
#############################################

######################################
###   Section 4.1 - UMAP & t-SNE   ###
######################################

library(Biobase)
library(Rtsne)
library(umap)

# update defaults for umap to contain a stable random_state (seed)
custom_umap <- umap::umap.defaults
custom_umap$random_state <- 42

# run UMAP
umap_out <-
  umap(t(log2(assayDataElement(target_data , elt = "q_norm"))),  
       config = custom_umap)

pData(target_data)[, c("UMAP1", "UMAP2")] <- umap_out$layout[, c(1,2)]

ggplot(pData(target_data),
       aes(x = UMAP1, y = UMAP2, color = region, shape = class)) +
  geom_point(size = 3, alpha = 0.5) +
  theme_bw()

ggsave(sprintf("%s/plots/4_1_unsupervised_analysis_umap.png", path), device='png')


# run tSNE
set.seed(42) # set the seed for tSNE as well

tsne_out <-
  Rtsne(t(log2(assayDataElement(target_data , elt = "q_norm"))),
        perplexity = ncol(target_data)*.15)

pData(target_data)[, c("tSNE1", "tSNE2")] <- tsne_out$Y[, c(1,2)]

ggplot(pData(target_data),
       aes(x = tSNE1, y = tSNE2, color = region, shape = class)) +
  geom_point(size = 3, alpha = 0.5) +
  theme_bw()

ggsave(sprintf("%s/plots/4_1_unsupervised_analysis_tsne.png", path), device='png')


##################################################
###   Section 4.2 - Clustering high CV Genes   ###
##################################################

library(pheatmap)

# create a log2 transform of the data for analysis
assayDataElement(object = target_data, elt = "log_q") <-
  assayDataApply(target_data, 2, FUN = log, base = 2, elt = "q_norm")

# create CV function
calc_CV <- function(x) {sd(x) / mean(x)}
CV_dat <- assayDataApply(target_data, elt = "log_q", MARGIN = 1, calc_CV)
# show the highest CD genes and their CV values
highest_cv <- data.frame(sort(CV_dat, decreasing = TRUE)[1:5])
highest_cv$Genes <- rownames(highest_cv)
rownames(highest_cv) <- NULL
names(highest_cv)[names(highest_cv) == 'sort.CV_dat..decreasing...TRUE..1.5.'] <- 'Coefficient of variation'
write.csv(highest_cv, sprintf("%s/data/4_2_table_highest_cd_genes.csv", path), row.names=FALSE)

# Identify genes in the top 3rd of the CV values
GOI <- names(CV_dat)[CV_dat > quantile(CV_dat, 0.8)]

png(sprintf("%s/plots/4_2_clustering_genes_coefficient_of_variation.png", path))
pheatmap(assayDataElement(target_data[GOI, ], elt = "log_q"),
         scale = "row", 
         show_rownames = FALSE, show_colnames = FALSE,
         border_color = NA,
         clustering_method = "average",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         breaks = seq(-3, 3, 0.05),
         color = colorRampPalette(c("purple3", "black", "yellow2"))(120),
         annotation_col = pData(target_data)[, c("class", "segment", "region")])
dev.off()

# Save image
save.image(sprintf('%s/image/4_unsupervised_analysis.RData', path))
