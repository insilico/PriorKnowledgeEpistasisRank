#
#
# Read and preprocess the Japan data from GEO site
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
Elist <- getGEO("GSE76826")

# Elicit gene expression
expr <- exprs(Elist[[1]])
dim(expr)
# Elicit the titles (MDD and HC)
titles <- pData(Elist[[1]])[, 1]
length(titles)
# assign the status labels
jap_class <- ifelse(grepl("HC", titles), 0, 1)
# Making sure that probIDs and rownames are identical
probID <- Elist$GSE76826_series_matrix.txt.gz@featureData@data$ID
identical(as.character(probID), rownames(expr))
# Substitute probID with gene symboles
gene.symbol <- Elist$GSE76826_series_matrix.txt.gz@featureData@data$GENE_SYMBOL
length(gene.symbol)
rownames(expr) <- gene.symbol
expr <- expr[-which(rownames(expr)==""), ]
dim(expr)

# z-transform in control samples
ctrl_idx <- which(jap_class==0)
ctrl_expr <- expr[, ctrl_idx]
z.trans_mat <- matrix(0,nrow = nrow(expr),ncol=ncol(expr))
for (i in 1:nrow(ctrl_expr)){
  z.trans_mat[i, ] <- pnorm(as.numeric(expr[i, ]),
                            mean(as.numeric(ctrl_expr[i, ])),
                            sd(as.numeric(ctrl_expr[i, ])))
}

# boxplot the raw and normalized data
par(mfrow=c(2, 1))
boxplot(as.data.frame(expr), range=0, main="Gunma gene expression",xlab="Subjects", ylab="Raw Expression")

boxplot(as.data.frame(z.trans_mat), range=0, main="",xlab="Subjects", ylab="pnorm Expression")


jap_expr <- t(z.trans_mat)
colnames(jap_expr) <- rownames(expr)
rownames(jap_expr) <- colnames(expr)

save(jap_expr, file = "jap_expr.RData")
save(jap_class, file = "jap_class.RData")

