#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
path = args[1]

load(sprintf('%s/image/6_unsupervised_analysis.RData', path))

###############################################
###   Section 7 - Differential Expression   ###
###############################################

###############################################
###   Section 7.1 - Within Slide Analysis   ###
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
###   Section 7.2 - Interpreting the results table   ###
########################################################

library(knitr)

table <- kable(subset(results, Gene %in% goi & Subset == "normal"), digits = 3,
      caption = "DE results for Genes of Interest",
      align = "lc", row.names = FALSE)

file_conn <- file(sprintf("%s/data/7_2_table_differential_expression_genes_of_interest_within_slide_analysis.txt", path))
writeLines(table, file_conn)
close(file_conn)


################################################
###   Section 7.3 - Between Slide Analysis   ###
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

kable(subset(results2, Gene %in% goi & Subset == "tubule"), digits = 3,
      caption = "DE results for Genes of Interest",
      align = "lc", row.names = FALSE)

file_conn <- file(sprintf("%s/data/7_3_table_differential_expression_genes_of_interest_between_slide_analysis.txt", path))
writeLines(table, file_conn)
close(file_conn)

# Save image
save.image(sprintf('%s/image/7_differential_expression.RData', path))

