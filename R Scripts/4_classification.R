#
#
# Classification of MDD datasets using nestedCV feature selection based on PageRank, Katz, 
# Epistasis-Katz and EpistasisRank centrality methods and xgboost classifier with and without prior knowledge.
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
library(xgboost)

setwd(".../path to directory")
source("helper_functions.R")
# load MDD microarray data
load("camb_fltr.expr.RData")
load("camb_sub_class.RData")
load("jap_fltr.expr.RData")
load("jap_class.RData")

# load reGAIN and IMP query network
load("camb_reGAIN.RData")
load("jap_reGAIN.RData")
load("impNetwork.RData")

######################################
#------ Without Prior Knowledge -----#
######################################
#----------- nestedCV feature selection-PageRank ----------#
page_woPK_features <- NULL
page_woPK_train_acc <- NULL
page_woPK_test_acc <- NULL
# create outer folds
outer_folds <- caret::createFolds(camb_class, 5, list = FALSE)
# create inner folds
for (i in 1:5){
  page_woPK_inner.features <- NULL
  inner_folds <- caret::createFolds(camb_class[outer_folds!=i], 5, list = TRUE)
  for (j in 1:length(inner_folds)){
    fold_idx <- which(outer_folds != i)[-inner_folds[[j]]] 
    # correlation matrix
    camb_expre.cor <- cor(camb_fltr.expr[fold_idx,])
    # co-expression to binary
    camb_Adj_partial <- ifelse(abs(camb_expre.cor) > .15, 1, 0)
    diag(camb_Adj_partial) <- 0
    # Adjacency matrix to graph
    camb_Adj_g <- graph.adjacency(camb_Adj_partial)
    # page_woPKRank centrality
    camb_page_woPK <- page.rank(camb_Adj_g)$vector
    camb_page_woPK.sort <- sort(camb_page_woPK, T)
    top_camb_page_woPK <- camb_page_woPK.sort[1:200]
    page_woPK_inner.features <- c(page_woPK_inner.features, top_camb_page_woPK)
    page_woPK_features <- c(page_woPK_features, top_camb_page_woPK)
    
  }
  camb_page_woPK.expr <- camb_fltr.expr[, unique(as.character(names(page_woPK_inner.features)))] 
  # Training
  camb_dtrain <- xgb.DMatrix(as.matrix(camb_page_woPK.expr), label = camb_class)
  camb_page_woPK_bst.model <- xgb.train(data=camb_dtrain, max.depth=1, eta = .01, nthread = 2, 
                                        nround=1, gamma = 0, min_child_weight = 1, max_delta_step = 1, 
                                        subsample = 0.5, objective = "binary:logistic")
  
  camb_page_woPK_bst.pred <- predict(camb_page_woPK_bst.model, as.matrix(camb_page_woPK.expr))
  page_woPK_train_acc <- c(page_woPK_train_acc, accuracy(camb_class, camb_page_woPK_bst.pred)$prop.correct)
  
  # Testing
  jap_page_woPK.expr <- jap_fltr.expr[, unique(as.character(names(page_woPK_inner.features)))]
  jap_page_woPK_bst.pred <- predict(camb_page_woPK_bst.model, as.matrix(jap_page_woPK.expr))
  page_woPK_test_acc <- c(page_woPK_test_acc, accuracy(jap_class, jap_page_woPK_bst.pred)$prop.correct)
}
save(page_woPK_features, file = "page_woPK_features.RData")
pr_wo_train <- page_woPK_train_acc[which.min(page_woPK_train_acc - page_woPK_test_acc)]
pr_wo_test <- page_woPK_test_acc[which.min(page_woPK_train_acc - page_woPK_test_acc)]

#------------- nestedCV feature selection-Katz -------------#
katz_woPK_features <- NULL
katz_woPK_train_acc <- NULL
katz_woPK_test_acc <- NULL
# create outer folds
outer_folds <- caret::createFolds(camb_class, 5, list = FALSE)
# create inner folds
for (i in 1:5){
  katz_woPK_inner.features <- NULL
  inner_folds <- caret::createFolds(camb_class[outer_folds!=i], 5, list = TRUE)
  for (j in 1:length(inner_folds)){
    fold_idx <- which(outer_folds != i)[-inner_folds[[j]]] 
    # correlation matrix
    camb_expre.cor <- cor(camb_fltr.expr[fold_idx,])
    # co-expression to binary
    camb_Adj_partial <- ifelse(abs(camb_expre.cor) > .15, 1, 0)
    diag(camb_Adj_partial) <- 0
    # katz_woPK centrality
    a <- eigen(camb_Adj_partial)
    alpha = signif(1/max(a$values), digits = 2)
    beta <- rep(1, nrow(camb_Adj_partial))/nrow(camb_Adj_partial)
    camb_katz_woPK <- katz.centrality(camb_Adj_partial, alpha, beta)
    names(camb_katz_woPK) <- colnames(camb_Adj_partial)
    camb_katz_woPK.sort <- sort(camb_katz_woPK, T)
    top_camb_katz_woPK <- camb_katz_woPK.sort[1:200]
    katz_woPK_inner.features <- c(katz_woPK_inner.features, top_camb_katz_woPK)
    katz_woPK_features <- c(katz_woPK_features, top_camb_katz_woPK)
  }
  camb_katz_woPK.expr <- camb_fltr.expr[, unique(as.character(names(katz_woPK_inner.features)))] 
  # Training
  camb_dtrain <- xgb.DMatrix(as.matrix(camb_katz_woPK.expr), label = camb_class)
  camb_katz_woPK_bst.model <- xgb.train(data=camb_dtrain, max.depth=1, eta = .01, nthread = 2, 
                                        nround=1, gamma = 0, min_child_weight = 1, max_delta_step = 1, 
                                        subsample = 0.5, objective = "binary:logistic")
  
  camb_katz_woPK_bst.pred <- predict(camb_katz_woPK_bst.model, as.matrix(camb_katz_woPK.expr))
  katz_woPK_train_acc <- c(katz_woPK_train_acc, accuracy(camb_class, camb_katz_woPK_bst.pred)$prop.correct)
  
  # Testing
  jap_katz_woPK.expr <- jap_fltr.expr[, unique(as.character(names(katz_woPK_inner.features)))]
  jap_katz_woPK_bst.pred <- predict(camb_katz_woPK_bst.model, as.matrix(jap_katz_woPK.expr))
  katz_woPK_test_acc <- c(katz_woPK_test_acc, accuracy(jap_class, jap_katz_woPK_bst.pred)$prop.correct)
}
save(katz_woPK_features, file = "katz_woPK_features.RData")
k_wo_train <- katz_woPK_train_acc[which.min(katz_woPK_train_acc - katz_woPK_test_acc)]
k_wo_test <- katz_woPK_test_acc[which.min(katz_woPK_train_acc - katz_woPK_test_acc)]

#---------- nestedCV feature selection-EpistasisKatz -----------#
ekatz_woPK_features <- NULL
ekatz_woPK_train_acc <- NULL
ekatz_woPK_test_acc <- NULL
# create outer folds
outer_folds <- caret::createFolds(camb_class, 5, list = FALSE)
# create inner folds
for (i in 1:5){
  ekatz_woPK_inner.features <- NULL
  inner_folds <- caret::createFolds(camb_class[outer_folds!=i], 5, list = TRUE)
  for (j in 1:length(inner_folds)){
    fold_idx <- which(outer_folds != i)[-inner_folds[[j]]]   
    # reGAIN
    camb_expr <- camb_fltr.expr[fold_idx,]
    class <- camb_class[fold_idx]
    regain_matrix <- data.frame(camb_expr, class)
    library(Rinbix)
    rownames(regain_matrix) <- NULL
    camb_Regain <- Rinbix::regainInbix(regain_matrix)
    # ekatz_woPK centrality
    alpha = 1/mean(colSums(camb_Regain$reGAIN))
    beta = diag(camb_Regain$reGAIN)
    diag(camb_Regain$reGAIN) <- 0
    camb_ekatz_woPK <- katz.centrality(camb_Regain$reGAIN, alpha, beta)
    names(camb_ekatz_woPK) <- colnames(camb_Regain$reGAIN)
    camb_ekatz_woPK.sort <- sort(camb_ekatz_woPK, T)
    top_camb_ekatz_woPK <- camb_ekatz_woPK.sort[1:200]
    ekatz_woPK_inner.features <- c(ekatz_woPK_inner.features, top_camb_ekatz_woPK)
    ekatz_woPK_features <- c(ekatz_woPK_features, top_camb_katz)
  }
  camb_ekatz_woPK.expr <- camb_fltr.expr[, unique(as.character(names(ekatz_woPK_inner.features)))] 
  # Training
  camb_dtrain <- xgb.DMatrix(as.matrix(camb_ekatz_woPK.expr), label = camb_class)
  camb_ekatz_woPK_bst.model <- xgb.train(data=camb_dtrain, max.depth=2, eta = .1, nthread = 2, 
                                         nround=1, gamma = 0, min_child_weight = 1, max_delta_step = 1, 
                                         subsample = 0.5, objective = "binary:logistic")
  
  camb_ekatz_woPK_bst.pred <- predict(camb_ekatz_woPK_bst.model, as.matrix(camb_ekatz_woPK.expr))
  ekatz_woPK_train_acc <- c(ekatz_woPK_train_acc, accuracy(camb_class, camb_ekatz_woPK_bst.pred)$prop.correct)
  
  # Testing
  jap_ekatz_woPK.expr <- jap_fltr.expr[, unique(as.character(names(ekatz_woPK_inner.features)))]
  jap_ekatz_woPK_bst.pred <- predict(camb_ekatz_woPK_bst.model, as.matrix(jap_ekatz_woPK.expr))
  ekatz_woPK_test_acc <- c(ekatz_woPK_test_acc, accuracy(jap_class, jap_ekatz_woPK_bst.pred)$prop.correct)
}
save(ekatz_woPK_features, file = "ekatz_woPK_features.RData")
ek_wo_train <- ekatz_woPK_train_acc[which.min(ekatz_woPK_train_acc - ekatz_woPK_test_acc)]
ek_wo_test <- ekatz_woPK_test_acc[which.min(ekatz_woPK_train_acc - ekatz_woPK_test_acc)]


#------------ nestedCV feature selection-EpistasisRank  ------------#
ER_woPK_features <- NULL
ER_woPK_train_acc <- NULL
ER_woPK_test_acc <- NULL
# create outer folds
outer_folds <- caret::createFolds(camb_class, 5, list = FALSE)
# create inner folds
for (i in 1:5){
  ER_woPK_inner.features <- NULL
  inner_folds <- caret::createFolds(camb_class[outer_folds!=i], 5, list = TRUE)
  for (j in 1:length(inner_folds)){
    fold_idx <- which(outer_folds != i)[-inner_folds[[j]]] 
    # reGAIN
    camb_expr <- camb_fltr.expr[fold_idx,]
    class <- camb_class[fold_idx]
    regain_matrix <- data.frame(camb_expr, class)
    library(Rinbix)
    rownames(regain_matrix) <- NULL
    camb_Regain <- regainInbix(regain_matrix)
    # ER_woPKrank centrality
    camb_ER_woPKrank <- snprankInbix(camb_Regain$reGAIN, gamma = .85)
    top_camb_ER_woPK <- camb_ER_woPKrank[1:200, ]
    ER_woPK_inner.features <- c(ER_woPK_inner.features, as.character(top_camb_ER_woPK$gene))
    ER_woPK_features <- c(ER_woPK_features, as.character(top_camb_ER_woPK$gene))
  }
  camb_ER_woPK.expr <- camb_fltr.expr[, unique(ER_woPK_inner.features)] 
  # Training
  camb_dtrain <- xgb.DMatrix(as.matrix(camb_ER_woPK.expr), label = camb_class)
  camb_ER_woPK_bst.model <- xgb.train(data=camb_dtrain, max.depth=1, eta = .01, nthread = 2, 
                                      nround=1, gamma = 0, min_child_weight = 1, max_delta_step = 1, 
                                      subsample = 0.5, objective = "binary:logistic")
  
  camb_ER_woPK_bst.pred <- predict(camb_ER_woPK_bst.model, as.matrix(camb_ER_woPK.expr))
  ER_woPK_train_acc <- c(ER_woPK_train_acc, accuracy(camb_class, camb_ER_woPK_bst.pred)$prop.correct)
  
  # Testing
  jap_ER_woPK.expr <- jap_fltr.expr[, unique(ER_woPK_inner.features)]
  jap_ER_woPK_bst.pred <- predict(camb_ER_woPK_bst.model, as.matrix(jap_ER_woPK.expr))
  ER_woPK_test_acc <- c(ER_woPK_test_acc, accuracy(jap_class, jap_ER_woPK_bst.pred)$prop.correct)
}
save(ER_woPK_features, file = "ER_woPK_features.RData")
er_wo_train <- ER_woPK_train_acc[which.min(ER_woPK_train_acc - ER_woPK_test_acc)]
er_wo_test <- ER_woPK_test_acc[which.min(ER_woPK_train_acc - ER_woPK_test_acc)]


#--------- nestedCV feature selection-coefficient of variantion ---------#
cov_features <- NULL
cov_train_acc <- NULL
cov_test_acc <- NULL
# create outer folds
outer_folds <- caret::createFolds(camb_class, 5, list = FALSE)
# create inner folds
for (i in 1:5){
  cov_inner.features <- NULL
  inner_folds <- caret::createFolds(camb_class[outer_folds!=i], 5, list = TRUE)
  for (j in 1:length(inner_folds)){
    fold_idx <- which(outer_folds != i)[-inner_folds[[j]]] 
    # CV filtering
    exprData <- t(camb_fltr.expr[fold_idx, ])
    thresh <- .5
    camb_expr.fltr_cov <- camb_fltr.expr
    while (ncol(camb_expr.fltr_cov)>200){
      camb_expr.fltr_cov <- t(cov.filter(exprData,threshold = thresh)$fdata)
      thresh <- thresh-0.00001
    }
    dim(camb_expr.fltr_cov)
    cov_inner.features <- c(cov_inner.features, colnames(camb_expr.fltr_cov))
    cov_features <- c(cov_features, colnames(camb_expr.fltr_cov))
  }
  camb_cov.expr <- camb_fltr.expr[, unique(as.character(cov_inner.features))] 
  # Training
  camb_dtrain <- xgb.DMatrix(as.matrix(camb_cov.expr), label = camb_class)
  camb_cov_bst.model <- xgb.train(data=camb_dtrain, max.depth=1, eta = .01, nthread = 2, 
                                  nround=1, gamma = 0, min_child_weight = 1, max_delta_step = 1, 
                                  subsample = 0.5, objective = "binary:logistic")
  
  camb_cov_bst.pred <- predict(camb_cov_bst.model, as.matrix(camb_cov.expr))
  cov_train_acc <- c(cov_train_acc, accuracy(camb_class, camb_cov_bst.pred)$prop.correct)
  
  # Testing
  jap_cov.expr <- jap_fltr.expr[, unique(as.character(cov_inner.features))]
  jap_cov_bst.pred <- predict(camb_cov_bst.model, as.matrix(jap_cov.expr))
  cov_test_acc <- c(cov_test_acc, accuracy(jap_class, jap_cov_bst.pred)$prop.correct)
}
save(cov_features, file = "cov_features.RData")
cov_train <- cov_train_acc[which.min(cov_train_acc - cov_test_acc)]
cov_test <- cov_test_acc[which.min(cov_train_acc - cov_test_acc)]


######################################
#------- With Prior Knowledge -------#
######################################
# prior knowledge
impgraph <- graph_from_adjacency_matrix(impNetwork, mode = "undirected", diag = TRUE, add.rownames = NA)
IMP_deg <- degree(impgraph, normalized = TRUE, loops = FALSE)

#----------- nestedCV feature selection-PageRank -----------#
page_wPK_features <- NULL
page_wPK_train_acc <- NULL
page_wPK_test_acc <- NULL
# create outer folds
outer_folds <- caret::createFolds(camb_class, 5, list = FALSE)
# create inner folds
for (i in 1:5){
  page_wPK_inner.features <- NULL
  inner_folds <- caret::createFolds(camb_class[outer_folds!=i], 5, list = TRUE)
  for (j in 1:length(inner_folds)){
    fold_idx <- which(outer_folds != i)[-inner_folds[[j]]] 
    # correlation matrix
    camb_expre.cor <- cor(camb_fltr.expr[fold_idx,])
    # co-expression to binary
    camb_Adj_partial <- ifelse(abs(camb_expre.cor) > .15, 1, 0)
    diag(camb_Adj_partial) <- 0
    # Adjacency matrix to graph
    camb_Adj_g <- graph.adjacency(camb_Adj_partial)
    # page_wPKRank centrality
    camb_page_wPK <- as.numeric(my.pagerank(camb_Adj_partial, IMP_deg))
    names(camb_page_wPK) <- rownames(camb_Adj_partial) 
    camb_page_wPK.sort <- sort(camb_page_wPK, T)
    top_camb_page_wPK <- camb_page_wPK.sort[1:200]
    page_wPK_inner.features <- c(page_wPK_inner.features, as.character(names(top_camb_page_wPK)))
    page_wPK_features <- c(page_wPK_features, as.character(names(top_camb_page_wPK)))
  }
  camb_page_wPK.expr <- camb_fltr.expr[, unique(page_wPK_inner.features)] 
  # Training
  camb_dtrain <- xgb.DMatrix(as.matrix(camb_page_wPK.expr), label = camb_class)
  camb_page_wPK_bst.model <- xgb.train(data=camb_dtrain, max.depth=1, eta = .01, nthread = 2, 
                                       nround=1, gamma = 0, min_child_weight = 1, max_delta_step = 1, 
                                       subsample = 0.5, objective = "binary:logistic")
  
  camb_page_wPK_bst.pred <- predict(camb_page_wPK_bst.model, as.matrix(camb_page_wPK.expr))
  page_wPK_train_acc <- c(page_wPK_train_acc, accuracy(camb_class, camb_page_wPK_bst.pred)$prop.correct)
  
  # Testing
  jap_page_wPK.expr <- jap_fltr.expr[, unique(page_wPK_inner.features)]
  jap_page_wPK_bst.pred <- predict(camb_page_wPK_bst.model, as.matrix(jap_page_wPK.expr))
  page_wPK_test_acc <- c(page_wPK_test_acc, accuracy(jap_class, jap_page_wPK_bst.pred)$prop.correct)
}
save(page_wPK_features, file = "page_wPK_features.RData")
pr_w_train <- page_wPK_train_acc[which.min(page_wPK_train_acc - page_wPK_test_acc)]
pr_w_test <- page_wPK_test_acc[which.min(page_wPK_train_acc - page_wPK_test_acc)]


#---------- nestedCV feature selection-Katz -----------#
katz_wPK_features <- NULL
katz_wPK_train_acc <- NULL
katz_wPK_test_acc <- NULL
# create outer folds
outer_folds <- caret::createFolds(camb_class, 5, list = FALSE)
# create inner folds
for (i in 1:5){
  katz_wPK_inner.features <- NULL
  inner_folds <- caret::createFolds(camb_class[outer_folds!=i], 5, list = TRUE)
  for (j in 1:length(inner_folds)){
    fold_idx <- which(outer_folds != i)[-inner_folds[[j]]] 
    # correlation matrix
    camb_expre.cor <- cor(camb_fltr.expr[fold_idx,])
    # co-expression to binary
    camb_Adj_partial <- ifelse(abs(camb_expre.cor) > .15, 1, 0)
    diag(camb_Adj_partial) <- 0
    # katz_wPK centrality
    a <- eigen(camb_Adj_partial)
    alpha <- signif(1/max(a$values), digits = 2)
    camb_katz_wPK <- katz.centrality(camb_Adj_partial, alpha, beta = IMP_deg)
    names(camb_katz_wPK) <- rownames(camb_Adj_partial) 
    camb_katz_wPK.sort <- sort(camb_katz_wPK, T)
    top_camb_katz_wPK <- camb_katz_wPK.sort[1:200]
    katz_wPK_inner.features <- c(katz_wPK_inner.features, as.character(names(top_camb_katz_wPK)))
    katz_wPK_features <- c(katz_wPK_features, as.character(names(top_camb_katz_wPK)))
  }
  camb_katz_wPK.expr <- camb_fltr.expr[, unique(katz_wPK_inner.features)] 
  # Training
  camb_dtrain <- xgb.DMatrix(as.matrix(camb_katz_wPK.expr), label = camb_class)
  camb_katz_wPK_bst.model <- xgb.train(data=camb_dtrain, max.depth=1, eta = .01, nthread = 2, 
                                       nround=1, gamma = 0, min_child_weight = 1, max_delta_step = 1, 
                                       subsample = 0.5, objective = "binary:logistic")
  
  camb_katz_wPK_bst.pred <- predict(camb_katz_wPK_bst.model, as.matrix(camb_katz_wPK.expr))
  katz_wPK_train_acc <- c(katz_wPK_train_acc, accuracy(camb_class, camb_katz_wPK_bst.pred)$prop.correct)
  
  # Testing
  jap_katz_wPK.expr <- jap_fltr.expr[, unique(katz_wPK_inner.features)]
  jap_katz_wPK_bst.pred <- predict(camb_katz_wPK_bst.model, as.matrix(jap_katz_wPK.expr))
  katz_wPK_test_acc <- c(katz_wPK_test_acc, accuracy(jap_class, jap_katz_wPK_bst.pred)$prop.correct)
}
save(katz_wPK_features, file = "katz_wPK_features.RData")
k_w_train <- katz_wPK_train_acc[which.min(katz_wPK_train_acc - katz_wPK_test_acc)]
k_w_test <- katz_wPK_test_acc[which.min(katz_wPK_train_acc - katz_wPK_test_acc)]


#--------- nestedCV feature selection-EpistasisKatz ---------#
ekatz_wPK_features <- NULL
ekatz_wPK_train_acc <- NULL
ekatz_wPK_test_acc <- NULL
acc_temp <- NULL
# create outer folds
outer_folds <- caret::createFolds(camb_class, 5, list = FALSE)
# create inner folds
for (i in 1:5){
  ekatz_wPK_inner.features <- NULL
  inner_folds <- caret::createFolds(camb_class[outer_folds!=i], 5, list = TRUE)
  for (j in 1:length(inner_folds)){
    fold_idx <- which(outer_folds != i)[-inner_folds[[j]]] 
    # reGAIN
    camb_expr <- camb_fltr.expr[fold_idx,]
    class <- camb_class[fold_idx]
    regain_matrix <- data.frame(camb_expr, class)
    library(Rinbix)
    rownames(regain_matrix) <- NULL
    camb_Regain <- regainInbix(regain_matrix)
    # ekatz_wPK centrality
    beta <- diag(camb_Regain$reGAIN)
    diag(camb_Regain$reGAIN) <- 0
    camb_ekatz_wPK <- katz.centrality(camb_Regain$reGAIN, IMP_deg, beta)
    names(camb_ekatz_wPK) <- colnames(camb_expr) 
    camb_ekatz_wPK.sort <- sort(camb_ekatz_wPK, T)
    top_camb_ekatz_wPK <- camb_ekatz_wPK.sort[1:200]
    ekatz_wPK_inner.features <- c(ekatz_wPK_inner.features, as.character(names(top_camb_ekatz_wPK)))
    ekatz_wPK_features <- c(ekatz_wPK_features, as.character(names(top_camb_ekatz_wPK)))
  }
  camb_ekatz_wPK.expr <- camb_fltr.expr[, unique(ekatz_wPK_inner.features)] 
  # Training
  camb_dtrain <- xgb.DMatrix(as.matrix(camb_ekatz_wPK.expr), label = camb_class)
  camb_ekatz_wPK_bst.model <- xgb.train(data=camb_dtrain, max.depth=2, eta = .1, nthread = 2, 
                                        nround=1, gamma = 0, min_child_weight = 1, max_delta_step = 1, 
                                        subsample = 0.5, objective = "binary:logistic")
  
  camb_ekatz_wPK_bst.pred <- predict(camb_ekatz_wPK_bst.model, as.matrix(camb_ekatz_wPK.expr))
  ekatz_wPK_train_acc <- c(ekatz_wPK_train_acc, accuracy(camb_class, camb_ekatz_wPK_bst.pred)$prop.correct)
  
  # Testing
  jap_ekatz_wPK.expr <- jap_fltr.expr[, unique(ekatz_wPK_inner.features)]
  jap_ekatz_wPK_bst.pred <- predict(camb_ekatz_wPK_bst.model, as.matrix(jap_ekatz_wPK.expr))
  ekatz_wPK_test_acc <- c(ekatz_wPK_test_acc, accuracy(jap_class, jap_ekatz_wPK_bst.pred)$prop.correct)
}
save(ekatz_wPK_features, file = "ekatz_wPK_features.RData")
ek_w_train <- ekatz_wPK_train_acc[which.min(ekatz_wPK_train_acc - ekatz_wPK_test_acc)]
ek_w_test <- ekatz_wPK_test_acc[which.min(ekatz_wPK_train_acc - ekatz_wPK_test_acc)]


#---------- nestedCV feature selection-EpistasisRank ---------#
ER_wPK_features <- NULL
ER_wPK_train_acc <- NULL
ER_wPK_test_acc <- NULL
# create outer folds
outer_folds <- caret::createFolds(camb_class, 5, list = FALSE)
# create inner folds
for (i in 1:5){
  ER_wPK_inner.features <- NULL
  inner_folds <- caret::createFolds(camb_class[outer_folds!=i], 5, list = TRUE)
  for (j in 1:length(inner_folds)){
    fold_idx <- which(outer_folds != i)[-inner_folds[[j]]] 
    # reGAIN
    camb_expr <- camb_fltr.expr[fold_idx,]
    class <- camb_class[fold_idx]
    regain_matrix <- data.frame(camb_expr, class)
    library(Rinbix)
    rownames(regain_matrix) <- NULL
    camb_Regain <- regainInbix(regain_matrix)
    # ER_wPKrank centrality
    camb_ER_wPKrank <- EpistasisRank(camb_Regain$reGAIN, IMP_deg)
    top_camb_ER_wPK <- camb_ER_wPKrank[1:200,]
    ER_wPK_inner.features <- c(ER_wPK_inner.features, as.character(top_camb_ER_wPK$gene))
    ER_wPK_features <- c(ER_wPK_features, as.character(top_camb_ER_wPK$gene))
  }
  camb_ER_wPK.expr <- camb_fltr.expr[, unique(ER_wPK_inner.features)] 
  # Training
  camb_dtrain <- xgb.DMatrix(as.matrix(camb_ER_wPK.expr), label = camb_class)
  camb_ER_wPK_bst.model <- xgb.train(data=camb_dtrain, max.depth=1, eta = .01, nthread = 2, 
                                     nround=1, gamma = 0, min_child_weight = 1, max_delta_step = 1, 
                                     subsample = 0.5, objective = "binary:logistic")
  
  camb_ER_wPK_bst.pred <- predict(camb_ER_wPK_bst.model, as.matrix(camb_ER_wPK.expr))
  ER_wPK_train_acc <- c(ER_wPK_train_acc, accuracy(camb_class, camb_ER_wPK_bst.pred)$prop.correct)
  
  # Testing
  jap_ER_wPK.expr <- jap_fltr.expr[, unique(ER_wPK_inner.features)]
  jap_ER_wPK_bst.pred <- predict(camb_ER_wPK_bst.model, as.matrix(jap_ER_wPK.expr))
  ER_wPK_test_acc <- c(ER_wPK_test_acc, accuracy(jap_class, jap_ER_wPK_bst.pred)$prop.correct)
}
save(ER_wPK_features, file = "ER_wPK_features.RData")
er_w_train <- ER_wPK_train_acc[which.min(ER_wPK_train_acc - ER_wPK_test_acc)]
er_w_test <- ER_wPK_test_acc[which.min(ER_wPK_train_acc - ER_wPK_test_acc)]

#--------- nestedCV feature selection-coefficient of variantion ---------#
cov_features <- NULL
cov_train_acc <- NULL
cov_test_acc <- NULL
# create outer folds
outer_folds <- caret::createFolds(camb_class, 5, list = FALSE)
# create inner folds
for (i in 1:5){
  cov_woPK_inner.features <- NULL
  inner_folds <- caret::createFolds(camb_class[outer_folds!=i], 5, list = TRUE)
  for (j in 1:length(inner_folds)){
    fold_idx <- which(outer_folds != i)[-inner_folds[[j]]] 
    # CV filtering
    exprData <- t(camb_fltr.expr[fold_idx, ])
    thresh <- .5
    camb_expr.fltr_cov <- camb_fltr.expr
    while (ncol(camb_expr.fltr_cov)>200){
      camb_expr.fltr_cov <- t(cov.filter(exprData,threshold = thresh)$fdata)
      thresh <- thresh-0.001
    }
    cov_inner.features <- c(cov_inner.features, colnames(camb_expr.fltr_cov))
    cov_features <- c(cov_features, colnames(camb_expr.fltr_cov))
  }
  camb_cov.expr <- camb_fltr.expr[, unique(as.character(cov_inner.features))] 
  # Training
  camb_dtrain <- xgb.DMatrix(as.matrix(camb_cov.expr), label = camb_class)
  camb_cov_bst.model <- xgb.train(data=camb_dtrain, max.depth=1, eta = .01, nthread = 2, 
                                  nround=1, gamma = 0, min_child_weight = 1, max_delta_step = 1, 
                                  subsample = 0.5, objective = "binary:logistic")
  
  camb_cov_bst.pred <- predict(camb_cov_bst.model, as.matrix(camb_cov.expr))
  cov_train_acc <- c(cov_train_acc, accuracy(camb_class, camb_cov_bst.pred)$prop.correct)
  
  # Testing
  jap_cov.expr <- jap_fltr.expr[, unique(as.character(cov_inner.features))]
  jap_cov_bst.pred <- predict(camb_cov_bst.model, as.matrix(jap_cov.expr))
  cov_test_acc <- c(cov_test_acc, accuracy(jap_class, jap_cov_bst.pred)$prop.correct)
}
save(cov_features, file = "cv_features.RData")
cov_train <- cov_train_acc[which.min(cov_train_acc - cov_test_acc)]
cov_test <- cov_test_acc[which.min(cov_train_acc - cov_test_acc)]


#--------- privateEC classificaton ---------#
camb_data <- data.frame(camb_fltr.expr, class = as.factor(camb_class))
jap_data <- data.frame(jap_fltr.expr, class = as.factor(jap_class))

train.splits <- splitDataset(camb_data, pct.train = 0.5, pct.holdout = 0.5, label = "class")
train.splits$train[, "class"] <- as.factor(ifelse(train.splits$train[, "class"] == -1, 1, 0))
train.splits$holdout[, "class"] <- as.factor(ifelse(train.splits$holdout[, "class"] == -1, 1, 0))
pec_results <- privateEC(train.ds = train.splits$train, holdout.ds = train.splits$holdout, validation.ds = jap_data,
                         label = "class", method.model = "classification", is.simulated = FALSE, tau.param = 100,
                         update.freq = 50, rf.ntree = 500, verbose = FALSE)

multiple.max.indices <- which(pec_results$algo.acc$holdout.acc==max(pec_results$algo.acc$holdout.acc))
last.max <- multiple.max.indices[length(multiple.max.indices)]
pec_train <- pec_results$algo.acc$holdout.acc[last.max]
pec_hold <- pec_results$algo.acc$validation.acc[last.max]


#--------- xgboost only classification ---------#
# Training
camb_dtrain <- xgb.DMatrix(as.matrix(camb_fltr.expr), label = camb_class)
camb_xgb_model <- xgb.train(data=camb_dtrain, max.depth=4, eta = .3, nthread = 2, 
                            nround=2, gamma = 0, min_child_weight = 1, max_delta_step = 1, 
                            subsample = 1, objective = "binary:logistic")

camb_xgb_pred <- predict(camb_xgb_model, as.matrix(camb_fltr.expr))
xgb_train <- accuracy(camb_class, camb_xgb_pred)$prop.correct

# Testing
jap_xgb_pred <- predict(camb_xgb_model, as.matrix(jap_fltr.expr))
xgb_test <- accuracy(jap_class, jap_xgb_pred)$prop.correct


###############################
#---- Classification plot ----#
###############################
PR.Katz_df <- data.frame(
  Data = factor(c("Cambridge-Train", "Japan-Test", "Cambridge-Train", "Japan-Test", 
                  "Cambridge-Train", "Japan-Test", "Cambridge-Train", "Japan-Test"), levels = c("Cambridge-Train", 
                                                                                                "Japan-Test")),
  Method = factor(c("PR", "PR", "Katz", "Katz", "PR", "PR", "Katz", "Katz"), 
                  levels = c("PR", "Katz")),
  Accuracy = c(pr_wo_train, pr_wo_test, k_wo_train, k_wo_test, pr_w_train, pr_w_test, k_w_train, k_w_test),
  Status  = c( "no Prior Knowledge", "no Prior Knowledge", "no Prior Knowledge", "no Prior Knowledge", 
               "Prior Knowledge", "Prior Knowledge", "Prior Knowledge", "Prior Knowledge")
)

EK.ER_df <- data.frame(
  Data = factor(c("Cambridge-Train", "Japan-Test", "Cambridge-Train", "Japan-Test", 
                  "Cambridge-Train", "Japan-Test", "Cambridge-Train", "Japan-Test"), levels = c("Cambridge-Train", 
                                                                                                "Japan-Test")),
  Method = factor(c("EK", "EK", "ER", "ER", 
                    "EK", "EK", "ER", "ER"), 
                  levels = c("EK", "ER")),
  Accuracy = c(ek_wo_train, ek_wo_test, er_wo_train, er_wo_test, ek_w_train, ek_w_test, er_w_train, er_w_test),
  Status = c("no Prior Knowledge", "no Prior Knowledge", "no Prior Knowledge", "no Prior Knowledge", 
             "Prior Knowledge", "Prior Knowledge", "Prior Knowledge", "Prior Knowledge")
)

cov_pec_df <- data.frame(
  Data = factor(c("Cambridge-Train", "Cambridge-Train", "Japan-Test", "Japan-Test")),
  Method = factor(c("CoV", "xgb", "CoV", "xgb"), levels = c("CoV", "xgb")),
  Accuracy = c(cov_train, xgb_train, cov_test, xgb_test),
  Status = c("no Prior Knowledge", "no Prior Knowledge", "no Prior Knowledge", "no Prior Knowledge")
)


p1 <- ggplot(PR.Katz_df, aes(x = Data, y = Accuracy, group = Method)) +
  ggtitle("Co-expression Networks") + theme(plot.title = element_text(hjust = 0.5)) +
  geom_point() + geom_line(aes(linetype=Method)) +
  geom_point(data = cov_pec_df, mapping = aes(x = Data, y = Accuracy, shape = Method)) +
  theme(axis.title.x=element_blank()) +
  scale_shape_manual(values=c(6,5)) +
  scale_linetype_manual(values=c("twodash", "solid"))+
  facet_wrap(~Status) +
  theme(legend.position=c(0.95, 0.5), legend.text=element_text(size=6), legend.title = element_text(size=6)) +
  guides(group = guide_legend(order=1),
         shape = guide_legend(order=2))
p2 <- ggplot(EK.ER_df, aes(x = Data, y = Accuracy, group = Method)) +
  ggtitle("Epistasis Networks") + theme(plot.title = element_text(hjust = 0.5)) +
  geom_point() + geom_line(aes(linetype=Method)) +
  geom_point(data = cov_pec_df, mapping = aes(x = Data, y = Accuracy, shape = Method)) +
  theme(axis.title.x=element_blank()) +
  scale_shape_manual(values=c(6,5)) +
  scale_linetype_manual(values=c("twodash", "solid"))+
  facet_wrap(~Status) +
  theme(legend.position=c(0.95, 0.5), legend.text=element_text(size=6), legend.title = element_text(size=6)) +
  guides(group = guide_legend(order=1),
         shape = guide_legend(order=2))
grid.arrange(p1, p2, nrow = 2)