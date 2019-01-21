#
# Pathway analysis
#
# Developed by: Saeid Parvandeh
#
# --------------------------------------------------------------
rm(list=ls())

library(igraph)
library(SDMTools)
library(Rinbix)
library(ggplot2)
library(caret)

#setwd(".../path to directory")
source("helper_functions.R")
# load MDD microarray data
load("camb_fltr.expr.RData")
load("camb_class_equalized.RData")
load("jap_fltr.expr.RData")
load("jap_class.RData")

# load reGAIN and IMP query network
load("camb_reGAIN.RData")
load("jap_reGAIN.RData")
load("impNetwork.RData")

# creating prior knowledge network
impgraph <- graph_from_adjacency_matrix(impNetwork, mode = "undirected", diag = TRUE, add.rownames = NA)
IMP_deg <- degree(impgraph, normalized = TRUE, loops = FALSE)

######## Cambridge Data set ##########
# ----------------------------------------------------
# co-expre+PageRank w/o prior knowledge
# ----------------------------------------------------
# correlation matrix
camb_expre.cor <- cor(camb_fltr.expr)
# co-expression to binary
camb_Adj <- ifelse(abs(camb_expre.cor) > .15, 1, 0)
diag(camb_Adj) <- 0

# PageRank_woPK centrality
camb_page_woPK <- Rinbix::PageRank(camb_Adj)[, 1]
camb_page_woPK.sort <- sort(camb_page_woPK, T)
top_camb_page_woPK <- camb_page_woPK.sort[1:200]


# ----------------------------------------------------
# co-expre+Katz w/o prior knowledge
# ----------------------------------------------------
a <- eigen(camb_Adj)
beta <- rep(1, nrow(camb_Adj))/nrow(camb_Adj)
alpha <- 1/max(a$values) - (1/max(a$values))/100
camb_katz_woPK <- Rinbix::EpistasisKatz(camb_Adj, alpha, beta)
names(camb_katz_woPK) <- colnames(camb_Adj)
camb_katz_woPK.sort <- sort(camb_katz_woPK, T)
top_camb_katz_woPK <- camb_katz_woPK.sort[1:200]

# ----------------------------------------------------
# reGAIN+EpistasisKatz w/o prior knowledge
# ----------------------------------------------------
alpha = 1/mean(colSums(camb_reGAIN))
beta = diag(camb_reGAIN)
camb_regain.Katz <- camb_reGAIN
diag(camb_regain.Katz) <- 0
camb_EK_woPK <- Rinbix::EpistasisKatz(camb_regain.Katz, alpha, beta)
names(camb_EK_woPK) <- rownames(camb_reGAIN) 
camb_EK_woPK.sort <- sort(camb_EK_woPK, T)
top_camb_EK_woPK <- camb_EK_woPK.sort[1:200]

# --------------------------------------------------
# reGAIN+EpistasisRank w/o prior knowledge
# --------------------------------------------------
camb_ER_woPK <- Rinbix::EpistasisRank(camb_reGAIN, Gamma_vec = .85)
top_camb_ER_woPK <- camb_ER_woPK[1:200, ]

# --------------------------------------------------
# Using CoV filtering
# --------------------------------------------------
exprData <- t(camb_fltr.expr)
top_camb_cov <- t(Rinbix::geneLowCoefOfVarFilter(exprData,.4874)$fdata) # 200 genes (sub = , ori = .4874)
dim(top_camb_cov)

# ----------------------------------------------------
# co-expre+PageRank w/prior knowledge
# ----------------------------------------------------
# correlation matrix
camb_expre.cor <- cor(camb_fltr.expr)
# co-expression to binary
camb_Adj <- ifelse(abs(camb_expre.cor) > .15, 1, 0)
diag(camb_Adj) <- 0

# PageRank_wPK centrality
camb_page_wPK <- as.numeric(Rinbix::PageRank(camb_Adj, IMP_deg))
names(camb_page_wPK) <- rownames(camb_Adj) 
camb_page_wPK.sort <- sort(camb_page_wPK, T)
top_camb_page_wPK <- camb_page_wPK.sort[1:200]

# ----------------------------------------------------
# co-expre+Katz w/prior knowledge
# ----------------------------------------------------
a <- eigen(camb_Adj)
1/max(a$values)
alpha <- 1/max(a$values) - (1/max(a$values))/100
camb_katz_wPK <- Rinbix::EpistasisKatz(camb_Adj, alpha, beta = IMP_deg)
names(camb_katz_wPK) <- rownames(camb_Adj) 
camb_katz_wPK.sort <- sort(camb_katz_wPK, T)
top_camb_katz_wPK <- camb_katz_wPK.sort[1:200]

# ----------------------------------------------------
# reGAIN+EpistasisKatz w/prior knowledge
# ----------------------------------------------------
beta <- diag(camb_reGAIN)
camb_regain.Katz <- camb_reGAIN
diag(camb_regain.Katz) <- 0
camb_EK_wPK <- Rinbix::EpistasisKatz(camb_regain.Katz, IMP_deg, beta)
names(camb_EK_wPK) <- rownames(camb_Adj) 
camb_EK_wPK.sort <- sort(camb_EK_wPK, T)
top_camb_EK_wPK <- camb_EK_wPK.sort[1:200]

# --------------------------------------------------
# reGAIN+EpistasisRank w/prior knowledge
# --------------------------------------------------
camb_ER_wPK <- Rinbix::EpistasisRank(camb_reGAIN, IMP_deg)
top_camb_ER_wPK <- camb_ER_wPK[1:200, ]

####################################
#--------- Pathway plot -----------#
####################################
pathways <- read.csv("pathways.csv", header = TRUE)

my_pathways_woPK <- cbind(names(top_camb_page_woPK), names(top_camb_katz_woPK), 
                          names(top_camb_EK_woPK), as.character(top_camb_ER_woPK$gene), colnames(top_camb_cov))
colnames(my_pathways_woPK) <- c("Page_woPK", "Katz_woPK", "EK_woPK", "ER_woPK", "CoV")
save(my_pathways_woPK, file = "my_pathways_woPK.RData")
pathway_mat_woPK <- matrix(0, 5, 9)
for (i in 1:5){
  for (j in 1:9){
    pathway_mat_woPK[i, j] <- length(intersect(my_pathways_woPK[, i], as.character(pathways[, j])))
  }
}

my_pathways_wPK <- cbind(names(top_camb_page_wPK), names(top_camb_katz_wPK), 
                         names(top_camb_EK_wPK), as.character(top_camb_ER_wPK$gene))

colnames(my_pathways_wPK) <- c("Page_wPK", "Katz_wPK", "EK_wPK", "ER_wPK")
save(my_pathways_wPK, file = "my_pathways_wPK.RData")
pathway_mat_wPK <- matrix(0, 5, 9)
for (i in 1:4){
  for (j in 1:9){
    pathway_mat_wPK[i, j] <- length(intersect(my_pathways_wPK[, i], as.character(pathways[, j])))
  }
}

plot_df <- data.frame(PK_Status=c("woPK", "woPK", "woPK", "woPK", "woPK", "wPK", "wPK", "wPK", "wPK", "wPK"), 
                      Method=c("PR", "Katz", "EK", "ER", "CoV", "PR", "Katz", "EK", "ER", "CoV"), 
                      rbind(pathway_mat_woPK, pathway_mat_wPK))

p1 <- ggplot(data=plot_df, aes(x=Method, y=X1, fill=PK_Status)) +
  labs(title="INHIBITION OF VOLTAGE GATED \n CD CA2 CHANNELS", x="Ranking Method", y = "Number of Overlap Genes")+
  theme(axis.title = element_text(size = .5)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_brewer(palette="Paired") + theme_minimal() +
  scale_fill_manual(values=c("lightgray", "black")) +
  theme_classic()

p2 <- ggplot(data=plot_df, aes(x=Method, y=X2, fill=PK_Status)) +
  labs(title="NEURONAL SYSTEM", x="Ranking Method", y = "Number of Overlap Genes")+
  theme(axis.title = element_text(size = .5)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_brewer(palette="Paired") + theme_minimal() +
  scale_fill_manual(values=c("lightgray", "black")) +
  theme_classic()

p3 <- ggplot(data=plot_df, aes(x=Method, y=X3, fill=PK_Status)) +
  labs(title="POTASSIUM CHANNELS", x="Ranking Method", y = "Number of Overlap Genes")+
  theme(axis.title = element_text(size = .5)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_brewer(palette="Paired") + theme_minimal() +
  scale_fill_manual(values=c("lightgray", "black")) +
  theme_classic()

p4 <- ggplot(data=plot_df, aes(x=Method, y=X4, fill=PK_Status)) +
  labs(title="TRANSMISSION ACROSS CHEMICAL \n SYNAPSES", x="Ranking Method", y = "Number of Overlap Genes")+
  theme(axis.title = element_text(size = .5)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_brewer(palette="Paired") + theme_minimal() +
  scale_fill_manual(values=c("lightgray", "black")) +
  theme_classic()

p5 <- ggplot(data=plot_df, aes(x=Method, y=X5, fill=PK_Status)) +
  labs(title="GPCR DOWNSTREAM SIGNALING", x="Ranking Method", y = "Number of Overlap Genes")+
  theme(axis.title = element_text(size = .5)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_brewer(palette="Paired") + theme_minimal() +
  scale_fill_manual(values=c("lightgray", "black")) +
  theme_classic()

p6 <- ggplot(data=plot_df, aes(x=Method, y=X6, fill=PK_Status)) +
  labs(title="SLC MEDICATED TRANSMEMBRANCE \n TRANSPORT", x="Ranking Method", y = "Number of Overlap Genes")+
  theme(axis.title = element_text(size = .5)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_brewer(palette="Paired") + theme_minimal() +
  scale_fill_manual(values=c("lightgray", "black")) +
  theme_classic()

p7 <- ggplot(data=plot_df, aes(x=Method, y=X7, fill=PK_Status)) +
  labs(title="NEUROTRANSMITTER RECEPTOR \n BINDING", x="Ranking Method", y = "Number of Overlap Genes")+
  theme(axis.title = element_text(size = .5)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_brewer(palette="Paired") + theme_minimal() +
  scale_fill_manual(values=c("lightgray", "black")) +
  theme_classic()

p8 <- ggplot(data=plot_df, aes(x=Method, y=X8, fill=PK_Status)) +
  labs(title="RAS ACTIVATION UPON CA2 INFLUX \n THROUGH NMDA RECEPTOR", x="Ranking Method", y = "Number of Overlap Genes")+
  theme(axis.title = element_text(size = .5)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_brewer(palette="Paired") + theme_minimal() +
  scale_fill_manual(values=c("lightgray", "black")) +
  theme_classic()

p9 <- ggplot(data=plot_df, aes(x=Method, y=X9, fill=PK_Status)) +
  labs(title ="SIGNALLING BY GPCR", x="Ranking Method", y = "Number of Overlap Genes")+
  theme(plot.title = element_text(size = .2)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_brewer(palette="Paired") + theme_minimal() +
  scale_fill_manual(values=c("lightgray", "black")) +
  theme_classic()

gridExtra::grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, nrow = 3, ncol = 3)

######## Japan Data set ##########
# ----------------------------------------------------
# co-expre+PageRank w/o prior knowledge
# ----------------------------------------------------
jap_expre.cor <- cor(jap_fltr.expr)
# co-expression to binary
jap_Adj <- ifelse(abs(jap_expre.cor) > .15, 1, 0)
diag(jap_Adj) <- 0
# Adjacency matrix to graph
jap_Adj_g <- graph.adjacency(jap_Adj)

jap_page_woPK <- Rinbix::PageRank(jap_Adj_g)$vector
jap_page_woPK.sort <- sort(jap_page_woPK, T)
top_jap_page_woPK <- jap_page_woPK.sort[1:200]

# ----------------------------------------------------
# co-expre+Katz w/o prior knowledge
# ----------------------------------------------------
a <- eigen(jap_Adj)
alpha <- 1/max(a$values) - (1/max(a$values))/100
beta <- rep(1, nrow(jap_Adj))/nrow(jap_Adj)
jap_katz_woPK <- Rinbix::EpistasisKatz(jap_Adj, alpha, beta)
names(jap_katz_woPK) <- colnames(jap_Adj)
jap_katz_woPK.sort <- sort(jap_katz_woPK, T)
top_jap_katz_woPK <- jap_katz_woPK.sort[1:200]

# ----------------------------------------------------
# reGAIN+EpistasisKatz w/o prior knowledge
# ----------------------------------------------------
a = 1/mean(colSums(jap_reGAIN))
beta = diag(jap_reGAIN)
jap_regain.Katz <- jap_reGAIN
diag(jap_regain.Katz) <- 0
jap_EK_woPK <- Rinbix::EpistasisKatz(jap_regain.Katz, a, beta)
names(jap_EK_woPK) <- rownames(jap_reGAIN) 
jap_EK_woPK.sort <- sort(jap_EK_woPK, T)
top_jap_EK_woPK <- jap_EK_woPK.sort[1:200]

# --------------------------------------------------
# reGAIN+EpistasisRank w/o prior knowledge
# --------------------------------------------------
jap_ER_woPK <- Rinbix::EpistasisRank(jap_reGAIN, Gamma_vec = .85)
top_jap_ER_woPK <- jap_ER_woPK[1:200, ]

# --------------------------------------------------
# Using CoV filtering
# --------------------------------------------------
exprData <- t(jap_fltr.expr)
jap_expr.fltr_cov <- t(Rinbix::geneLowCoefOfVarFilter(exprData,.3795)$fdata) # 200 genes
dim(jap_expr.fltr_cov)

# ----------------------------------------------------
# co-expre+PageRank w/prior knowledge
# ----------------------------------------------------
jap_page_wPK <- as.numeric(Rinbix::PageRank(jap_Adj, IMP_deg))
names(jap_page_wPK) <- rownames(jap_Adj) 
jap_page_wPK.sort <- sort(jap_page_wPK, T)
top_jap_page_wPK <- jap_page_wPK.sort[1:200]

# ----------------------------------------------------
# co-expre+Katz w/prior knowledge
# ----------------------------------------------------
a <- eigen(jap_Adj)
alpha <- 1/max(a$values) - (1/max(a$values))/100
jap_katz_wPK <- Rinbix::EpistasisKatz(jap_Adj, alpha, beta = IMP_deg)
names(jap_katz_wPK) <- rownames(jap_Adj) 
jap_katz_wPK.sort <- sort(jap_katz_wPK, T)
top_jap_katz_wPK <- jap_katz_wPK.sort[1:200]

# ----------------------------------------------------
# reGAIN+EpistasisKatz w/prior knowledge
# ----------------------------------------------------
beta = diag(jap_reGAIN)
jap_regain.Katz <- jap_reGAIN
diag(jap_regain.Katz) <- 0
jap_EK_wPK <- Rinbix::EpistasisKatz(jap_regain.Katz, IMP_deg, beta)
names(jap_EK_wPK) <- rownames(jap_reGAIN) 
jap_EK_wPK.sort <- sort(jap_EK_wPK, T)
top_jap_EK_wPK <- jap_EK_wPK.sort[1:200]

# --------------------------------------------------
# reGAIN+EpistasisRank w/prior knowledge
# --------------------------------------------------
jap_ER_wPK <- Rinbix::EpistasisRank(jap_reGAIN, IMP_deg)
top_jap_ER_wPK <- jap_ER_wPK[1:200,]


####################################
#--------- Pathway plot -----------#
####################################
pathways <- read.csv("pathways.csv", header = TRUE)

my_pathways_woPK <- cbind(names(top_jap_page_woPK), names(top_jap_katz_woPK), 
                          names(top_jap_EK_woPK), as.character(top_jap_ER_woPK$gene), colnames(jap_expr.fltr_cov))
colnames(my_pathways_woPK) <- c("Page_woPK", "Katz_woPK", "EK_woPK", "ER_woPK", "CoV")
save(my_pathways_woPK, file = "my_pathways_woPK.RData")
pathway_mat_woPK <- matrix(0, 5, 9)
for (i in 1:5){
  for (j in 1:9){
    pathway_mat_woPK[i, j] <- length(intersect(my_pathways_woPK[, i], as.character(pathways[, j])))
  }
}

my_pathways_wPK <- cbind(names(top_jap_page_wPK), names(top_jap_katz_wPK), 
                         names(top_jap_EK_wPK), as.character(top_jap_ER_wPK$gene))

colnames(my_pathways_wPK) <- c("Page_wPK", "Katz_wPK", "EK_wPK", "ER_wPK")
save(my_pathways_wPK, file = "my_pathways_wPK.RData")
pathway_mat_wPK <- matrix(0, 5, 9)
for (i in 1:4){
  for (j in 1:9){
    pathway_mat_wPK[i, j] <- length(intersect(my_pathways_wPK[, i], as.character(pathways[, j])))
  }
}

plot_df <- data.frame(PK_Status=c("woPK", "woPK", "woPK", "woPK", "woPK", "wPK", "wPK", "wPK", "wPK", "wPK"), 
                      Method=c("PR", "Katz", "EK", "ER", "CoV", "PR", "Katz", "EK", "ER", "CoV"), 
                      rbind(pathway_mat_woPK, pathway_mat_wPK))

p1 <- ggplot(data=plot_df, aes(x=Method, y=X1, fill=PK_Status)) +
  labs(title="INHIBITION OF VOLTAGE GATED \n CD CA2 CHANNELS", x="Ranking Method", y = "Number of Overlap Genes")+
  theme(axis.title = element_text(size = .5)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_brewer(palette="Paired") + theme_minimal() +
  scale_fill_manual(values=c("lightgray", "black")) +
  theme_classic()

p2 <- ggplot(data=plot_df, aes(x=Method, y=X2, fill=PK_Status)) +
  labs(title="NEURONAL SYSTEM", x="Ranking Method", y = "Number of Overlap Genes")+
  theme(axis.title = element_text(size = .5)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_brewer(palette="Paired") + theme_minimal() +
  scale_fill_manual(values=c("lightgray", "black")) +
  theme_classic()

p3 <- ggplot(data=plot_df, aes(x=Method, y=X3, fill=PK_Status)) +
  labs(title="POTASSIUM CHANNELS", x="Ranking Method", y = "Number of Overlap Genes")+
  theme(axis.title = element_text(size = .5)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_brewer(palette="Paired") + theme_minimal() +
  scale_fill_manual(values=c("lightgray", "black")) +
  theme_classic()

p4 <- ggplot(data=plot_df, aes(x=Method, y=X4, fill=PK_Status)) +
  labs(title="TRANSMISSION ACROSS CHEMICAL \n SYNAPSES", x="Ranking Method", y = "Number of Overlap Genes")+
  theme(axis.title = element_text(size = .5)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_brewer(palette="Paired") + theme_minimal() +
  scale_fill_manual(values=c("lightgray", "black")) +
  theme_classic()

p5 <- ggplot(data=plot_df, aes(x=Method, y=X5, fill=PK_Status)) +
  labs(title="GPCR DOWNSTREAM SIGNALING", x="Ranking Method", y = "Number of Overlap Genes")+
  theme(axis.title = element_text(size = .5)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_brewer(palette="Paired") + theme_minimal() +
  scale_fill_manual(values=c("lightgray", "black")) +
  theme_classic()

p6 <- ggplot(data=plot_df, aes(x=Method, y=X6, fill=PK_Status)) +
  labs(title="SLC MEDICATED TRANSMEMBRANCE \n TRANSPORT", x="Ranking Method", y = "Number of Overlap Genes")+
  theme(axis.title = element_text(size = .5)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_brewer(palette="Paired") + theme_minimal() +
  scale_fill_manual(values=c("lightgray", "black")) +
  theme_classic()

p7 <- ggplot(data=plot_df, aes(x=Method, y=X7, fill=PK_Status)) +
  labs(title="NEUROTRANSMITTER RECEPTOR \n BINDING", x="Ranking Method", y = "Number of Overlap Genes")+
  theme(axis.title = element_text(size = .5)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_brewer(palette="Paired") + theme_minimal() +
  scale_fill_manual(values=c("lightgray", "black")) +
  theme_classic()

p8 <- ggplot(data=plot_df, aes(x=Method, y=X8, fill=PK_Status)) +
  labs(title="RAS ACTIVATION UPON CA2 INFLUX \n THROUGH NMDA RECEPTOR", x="Ranking Method", y = "Number of Overlap Genes")+
  theme(axis.title = element_text(size = .5)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_brewer(palette="Paired") + theme_minimal() +
  scale_fill_manual(values=c("lightgray", "black")) +
  theme_classic()

p9 <- ggplot(data=plot_df, aes(x=Method, y=X9, fill=PK_Status)) +
  labs(title ="SIGNALLING BY GPCR", x="Ranking Method", y = "Number of Overlap Genes")+
  theme(plot.title = element_text(size = .2)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_brewer(palette="Paired") + theme_minimal() +
  scale_fill_manual(values=c("lightgray", "black")) +
  theme_classic()

gridExtra::grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, nrow = 3, ncol = 3)



