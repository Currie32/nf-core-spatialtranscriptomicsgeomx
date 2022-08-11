#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
path = args[1]

load(sprintf('%s/image/4_unsupervised_analysis.RData', path))

###############################################
###   Section 5 - Differential Expression   ###
###############################################

###############################################
###   Section 5.1 - Within Slide Analysis   ###
###############################################

library(Biobase)

# convert test variables to factors
pData(target_data)$testRegion <- 
  factor(pData(target_data)$region, c("glomerulus", "tubule"))
pData(target_data)[["slide"]] <- 
  factor(pData(target_data)[["slide name"]])
assayDataElement(object = target_data, elt = "log_q") <-
  assayDataApply(target_data, 2, FUN = log, base = 2, elt = "q_norm")

# run LMM:
# formula follows conventions defined by the lme4 package
results <- c()
for(status in c("DKD", "normal")) {
  ind <- pData(target_data)$class == status
  mixedOutmc <-
    mixedModelDE(target_data[, ind],
                 elt = "log_q",
                 modelFormula = ~ testRegion + (1 + testRegion | slide),
                 groupVar = "testRegion",
                 nCores = parallel::detectCores(),
                 multiCore = FALSE)
  
  # format results as data.frame
  r_test <- do.call(rbind, mixedOutmc["lsmeans", ])
  tests <- rownames(r_test)
  r_test <- as.data.frame(r_test)
  r_test$Contrast <- tests
  
  # use lapply in case you have multiple levels of your test factor to
  # correctly associate gene name with it's row in the results table
  r_test$Gene <- 
    unlist(lapply(colnames(mixedOutmc),
                  rep, nrow(mixedOutmc["lsmeans", ][[1]])))
  r_test$Subset <- status
  r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
  r_test <- r_test[, c("Gene", "Subset", "Contrast", "Estimate", 
                       "Pr(>|t|)", "FDR")]
  results <- rbind(results, r_test)
}


########################################################
###   Section 5.2 - Interpreting the results table   ###
########################################################

library(knitr)
library(dplyr)
library(tidyr)

results_df <- data.frame(results)
results_df <- results_df[
  (results_df$Gene %in% goi & results_df$Subset == "normal"),
]
results_df <- results_df %>% mutate_if(is.numeric, round, digits=3)
names(results_df)[names(results_df) == 'Pr...t..'] <- 'P-value'

rownames(results_df) <- NULL

write.csv(
  results_df,
  sprintf("%s/data/5_2_table_differential_expression_genes_of_interest_within_slide_analysis.csv", path),
  row.names=FALSE
)


################################################
###   Section 5.3 - Between Slide Analysis   ###
################################################

# convert test variables to factors
pData(target_data)$testClass <-
  factor(pData(target_data)$class, c("normal", "DKD"))

# run LMM:
# formula follows conventions defined by the lme4 package
results2 <- c()
for(region in c("glomerulus", "tubule")) {
  ind <- pData(target_data)$region == region
  mixedOutmc <-
    mixedModelDE(target_data[, ind],
                 elt = "log_q",
                 modelFormula = ~ testClass + (1 | slide),
                 groupVar = "testClass",
                 nCores = parallel::detectCores(),
                 multiCore = FALSE)
  
  # format results as data.frame
  r_test <- do.call(rbind, mixedOutmc["lsmeans", ])
  tests <- rownames(r_test)
  r_test <- as.data.frame(r_test)
  r_test$Contrast <- tests
  
  # use lapply in case you have multiple levels of your test factor to
  # correctly associate gene name with it's row in the results table
  r_test$Gene <- 
    unlist(lapply(colnames(mixedOutmc),
                  rep, nrow(mixedOutmc["lsmeans", ][[1]])))
  r_test$Subset <- region
  r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
  r_test <- r_test[, c("Gene", "Subset", "Contrast", "Estimate", 
                       "Pr(>|t|)", "FDR")]
  results2 <- rbind(results2, r_test)
}

results2_df <- data.frame(results2)
results2_df <- results2_df[
  (results2_df$Gene %in% goi & results2_df$Subset == "tubule"),
]
results2_df <- results2_df %>% mutate_if(is.numeric, round, digits=3)
names(results2_df)[names(results2_df) == 'Pr...t..'] <- 'P-value'
rownames(results2_df) <- NULL

write.csv(
  results2_df,
  sprintf("%s/data/5_3_table_differential_expression_genes_of_interest_between_slide_analysis.csv", path),
  row.names=FALSE
)


# Save image
save.image(sprintf('%s/image/5_differential_expression.RData', path))
