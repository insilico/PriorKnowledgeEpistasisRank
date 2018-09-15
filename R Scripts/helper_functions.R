# Function: co-expression threshold estimation
thresh_fun <- function(thresh.points,wei.matrix)
{
  gAdjBool <- ifelse(wei.matrix > thresh.points, 1, 0)
  diag(gAdjBool) <- 0
  g.DB <- graph.adjacency(gAdjBool, "undirected")
  edges <- ecount(g.DB)
  full <- (ncol(gAdjBool)*(ncol(gAdjBool)-1))/2
  return (edges/full) 
}

# Function: centrality density estimation
density_fun <- function(thresh.points, wei.matrix){
  Adj <- ifelse(wei.matrix>thresh.points, 1, 0)
  diag(Adj) <- 0
  Adj_g <- graph.adjacency(Adj)
  deg <- degree(Adj_g)
  return (deg)
}

# Katz centrality
katz.centrality = function(A, alpha, beta) {
  n = nrow(A)
  I = diag(1, n); 
  c = solve(I - alpha * A,beta)  
  return(as.vector(c))
}

# Page rank centrality
my.pagerank <- function(A, damping){
  #make sure A is matrix type
  n <-nrow(A)
  rand.jump.col <- (1-damping)/n*matrix(rep(1, n),nrow=n)
  gDegs <- sapply(1:n, function(x) sum(A[x,]))
  Dk_cols<-matrix(rep(gDegs,n), nrow=n, byrow=TRUE)
  Amarkov <- A/Dk_cols
  my.sol <- solve(diag(n) - damping*Amarkov, rand.jump.col)
  return(my.sol)
}

# Modified SNPrank function
EpistasisRank <- function (G, Gamma_vec)
{
  n <- nrow(G)
  geneNames <- colnames(G)
  Gdiag <- diag(G)
  Gtrace <- sum(Gdiag)
  colsum <- colSums(G)
  diag(G) <- 0
  Gtrace <- Gtrace * n
  colsumG <- colSums(G)
  rowSumG <- rowSums(G)
  rowsum_denom <- matrix(0, n, 1)
  for (i in 1:n) {
    localSum <- 0
    for (j in 1:n) {
      factor <- ifelse(G[i, j] != 0, 1, 0)
      localSum <- localSum + (factor * colsumG[j])
    }
    rowsum_denom[i] <- localSum
  }
  if (length(Gamma_vec)==1){
    gamma_vec <- rep(Gamma_vec, n)
  } else {
    gamma_vec <- Gamma_vec
  }
  gamma_matrix <- matrix(nrow = n, ncol = n, data = rep(gamma_vec, n))
  if (Gtrace) {
    b <- ((1 - gamma_vec)/n) + (Gdiag/Gtrace)
  }
  else {
    b <- ((1 - gamma_vec)/n)
  }
  D <- matrix(nrow = n, ncol = n, data = c(0))
  diag(D) <- 1/colsumG
  I <- diag(n)
  temp <- I - gamma_matrix * G %*% D
  r <- solve(temp, b)
  snpranks <- r/sum(r)
  saveTable <- data.frame(gene = geneNames, snprank = snpranks)
  sortedTable <- saveTable[order(saveTable$snprank, decreasing = TRUE), ]
}

# Coefficient of Variation (C.V)
cov.filter <- function(dataMatrix, threshold) {
  # coefficient of variation filter
  mask <- apply(dataMatrix, 1, function(x) {(sd(x)/abs(mean(x))) < threshold})
  fdata <- dataMatrix[mask, ]
  # return the row mask and filtered data
  list(mask=mask, fdata=fdata)
}