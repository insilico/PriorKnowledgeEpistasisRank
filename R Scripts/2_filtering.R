#
#
# CoV filtering across datasets, Balance the samples, IMP query, reGAIN
#
# Author: Saeid Parvandeh
#
# -----------------------------------------------------------

rm(list=ls())

library(igraph)
library(SDMTools)
library(infotheo)
library(gdata)
library(MASS)

setwd(".../path to the directory")
source("helper_functions.R")
# load MDD microarray data
load("camb_expr.RData")
load("camb_class.RData")
load("jap_expr.RData")
load("jap_class.RData")

# balance the number of case/control subjects.
# ---------------------------------------------------
# case_idx <- which(camb_class==1)
# set.seed(123)
# half_case_idx <- sample(case_idx, length(case_idx)/2)
# camb_class <- camb_class[-half_case_idx]
# camb_expr <- camb_expr[-half_case_idx, ]
# save(camb_expr, file = "camb_expr_equalized.RData")
# save(camb_class, file = "camb_class_equalized.RData")
load("camb_expr_equalized.RData")
load("camb_class_equalized.RData")
# subsampling to the size of Japan data
dim(camb_expr)
table(camb_class)
# case
case_idx <- which(camb_class==1)
set.seed(456)
subsample_case_idx <- sample(case_idx, length(case_idx)-44)
# control
control_idx <- which(camb_class==0)
subsample_control_idx <- sample(control_idx, length(control_idx)-44)
camb_class <- camb_class[c(subsample_control_idx, subsample_case_idx)]
save(camb_class, file = "camb_sub_class.RData")
camb_expr <- camb_expr[c(subsample_control_idx,subsample_case_idx), ]
save(camb_expr, file = "camb_sub_expr.RData")

dim(camb_expr)
table(camb_class)
# ---------------------------------------------------

# Coefficient of Variation (CoV) filtering -- Cambridge, Japan
exprData1 <- t(camb_expr)
camb_expr.fltr_cov <- data.frame(t(Rinbix::geneLowCoefOfVarFilter(exprData1,.559095)$fdata)) 
dim(camb_expr.fltr_cov)

exprData2 <- t(jap_expr)
jap_expr.fltr_cov <- data.frame(t(Rinbix::geneLowCoefOfVarFilter(exprData2,.559095)$fdata)) 
dim(jap_expr.fltr_cov)

camb_genes <- as.character(colnames(camb_expr.fltr_cov))
jap_genes <- as.character(colnames(jap_expr.fltr_cov))

common_genes <- intersect(camb_genes, jap_genes)
length(common_genes)
camb_fltr.expr <- camb_expr.fltr_cov[, which(camb_genes %in% common_genes)]
camb_fltr.expr <- camb_fltr.expr[, order(colnames(camb_fltr.expr))]
dim(camb_fltr.expr)
save(camb_fltr.expr, file = "camb_fltr.expr.RData")
jap_fltr.expr <- jap_expr.fltr_cov[, which(jap_genes %in% common_genes)]
jap_fltr.expr <- jap_fltr.expr[, order(colnames(jap_fltr.expr))]
dim(jap_fltr.expr)
save(jap_fltr.expr, file = "jap_fltr.expr.RData")

# IMP query
# wirte all genes in a file
write.matrix(colnames(camb_fltr.expr), file = "topgenes.tab", sep = "\t")
# read the IMP query matrix
# please refer to readme file for computing IMP query matrix
impNetwork <- as.matrix(read.table("camb_jap_sub_impquery.adjmat", header = TRUE, sep = "\t"))
# replace all real values with 1, so that we will have a binary matrix
impNetwork[impNetwork > 0] <- 1
impNetwork <- matrix(as.integer(impNetwork), nrow = 5000)
rownames(impNetwork) <- colnames(camb_fltr.expr)
colnames(impNetwork) <- colnames(camb_fltr.expr)
save(impNetwork, file = "impNetwork.RData")

#### ---- Cambridge ---- ####
# inbix-reGAIN 
# prepare data to compute reGAIN
camb_regain_pheno <- data.frame(IID = rownames(camb_expr), FID = rownames(camb_expr), class = camb_class)
camb_regain_expr <- data.frame(IID = rownames(camb_expr), FID = rownames(camb_expr), camb_fltr.expr)
rownames(camb_regain_expr) <- NULL
write.matrix(camb_regain_pheno, file = "camb_pheno", sep = "\t")
write.matrix(camb_regain_expr, file = "camb_expre.num", sep = "\t")
# read reGAIN matrix
# please refer to readme file for computing reGAIN matrix
camb_reGAIN <- as.matrix(read.table("camb_sub.block.regain", header = T, sep = "\t"))
rownames(camb_reGAIN) <- colnames(camb_reGAIN)
save(camb_reGAIN, file = "camb_reGAIN.RData")
# alternative way to compute reGAIN
camb_data <- data.frame(camb_fltr.expr, class = camb_class)
camb_reGAIN <- regainInbix(camb_data)$reGAIN

#### ---- Japan ---- ####
# inbix-reGAIN 
# prepare data to compute reGAIN
jap_regain_pheno <- data.frame(IID = rownames(jap_expr), FID = rownames(jap_expr), class = jap_class)
jap_regain_expr <- data.frame(IID = rownames(jap_expr), FID = rownames(jap_expr), jap_fltr.expr)
rownames(jap_regain_expr) <- NULL
write.matrix(jap_regain_pheno, file = "jap_pheno", sep = "\t")
write.matrix(jap_regain_expr, file = "jap_expre.num", sep = "\t")
# read reGAIN matrix
# please refer to readme file for computing reGAIN matrix
jap_reGAIN <- as.matrix(read.table("jap.block.regain", header = T, sep = "\t"))
rownames(jap_reGAIN) <- colnames(jap_reGAIN)
save(jap_reGAIN, file = "jap_reGAIN.RData")
# alternative way to compute reGAIN
jap_data <- data.frame(jap_fltr.expr, class = jap_class)
jap_reGAIN <- regainInbix(jap_data)$reGAIN
