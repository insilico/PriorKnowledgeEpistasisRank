#
#
# Read and preprocess the Cambridge data from GEO site
#
# Author: Saeid Parvandeh
#
#----------------------------------------
rm(list=ls())
# read data using read.ilmn function
library(GEOquery)
library(Biobase)
# biocLite("preprocessCore")
library(preprocessCore)
# get data from GEO
Elist <- getGEO("GSE98793")

# Elicit gene expression
expr <- exprs(Elist[[1]])
dim(expr)
# Elicit the titles (MDD and CON)
titles <- pData(Elist[[1]])[, 1]
length(titles)
# assign the status labels
camb_class <- ifelse(grepl("control", titles), 0, 1)
# Making sure that probIDs and rownames are identical
probID <- Elist$GSE98793_series_matrix.txt.gz@featureData@data$ID
identical(as.character(probID), rownames(expr))
# Substitute probID with gene symboles
gene.symbol <- Elist$GSE98793_series_matrix.txt.gz@featureData@data$`Gene Symbol`
length(gene.symbol)
rownames(expr) <- gene.symbol
expr <- expr[-which(rownames(expr)==""), ]
dim(expr)

# z-transformation on control samples
ctrl_idx <- which(camb_class==0)
ctrl_expr <- expr[, ctrl_idx]
z.trans_mat <- matrix(0,nrow = nrow(expr),ncol=ncol(expr))
for (i in 1:nrow(ctrl_expr)){
  z.trans_mat[i, ] <- pnorm(as.numeric(expr[i, ]),
                            mean(as.numeric(ctrl_expr[i, ])),
                            sd(as.numeric(ctrl_expr[i, ])))
}

# boxplot the raw and normalized data
par(mfrow=c(2, 1))
boxplot(as.data.frame(expr), range=0, main="Cambridge gene expression",xlab="Subjects", ylab="Raw Expression")

boxplot(as.data.frame(z.trans_mat), range=0, main="",xlab="Subjects", ylab="pnorm Expression")


camb_expr <- t(z.trans_mat)
colnames(camb_expr) <- rownames(expr)
rownames(camb_expr) <- colnames(expr)

save(camb_expr, file = "camb_expr.RData")
save(camb_class, file = "camb_class.RData")

