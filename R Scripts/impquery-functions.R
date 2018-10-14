# impquery-functions.R - Bill White (bcw) - 8/30/17
#
# A port of Ahwan Pandey (ap)'s Python script impQuery.py - 11/11/12
# Functions to support impqery.R driver script.

library(RMariaDB)

# constants
INTERACTIONS_TABLE_GENE1_FIELD_IDX <- 1
INTERACTIONS_TABLE_GENE2_FIELD_IDX <- 2
INTERACTIONS_TABLE_CONFIDENCE_FIELD_IDX <- 3

RETRY_LIMIT <- 5

fillRegain <- function(impDbConn = NULL,
                       geneList = NULL,
                       scoreList = NULL,
                       hasMainEffectScores = FALSE,
                       ranked = FALSE,
                       confidenceThreshold = 0,
                       addMode = FALSE,
                       sifFilename = "",
                       verbose = FALSE) {
  # check function arguments
  if (is.null(impDbConn)) {
    stop("impDbConn is a required parameter")
  }
  if (is.null(geneList)) {
    stop("geneList is a required parameter")
  }
  if (length(geneList) < 1) {
    stop("geneList is empty")
  }
  inGeneListSize <- length(geneList)
  if (verbose) cat("-------------------------------------------------------------\n")
  if (verbose) cat("Begin fillRegain\n")
  if (verbose) cat("The gene list has", inGeneListSize, "\n")
  writeSifFile <- FALSE
  if (sifFilename == "") {
    if (verbose) cat("Write SIF file disabled\n")
    writeSifFile <- FALSE
  } else {
    if (verbose) cat("Write SIF file enabled\n")
    writeSifFile <- TRUE
  }
  sifFileCon <- NULL
  if (writeSifFile) {
    cat("Saving SIF file:", sifFilename, "\n")
    sifFileCon <- file(sifFilename, "w")
  }
  # fill a matrix with all combinations of pairs of gene indices
  regainMatrix <- matrix(nrow = inGeneListSize, ncol = inGeneListSize, data = c(0))
  combList <- combn(inGeneListSize, 2, simplify = FALSE)
  resultsList <- lapply(combList, function(x) {
    gene1Idx <- x[1]
    gene2Idx <- x[2]
    gene1Name <- geneList[gene1Idx]
    gene2Name <- geneList[gene2Idx]
    sqlString1 <- paste("select * from interactions where gene1='", gene1Name,
                        "' and gene2='", gene2Name, "'", sep = "")
    if (verbose) cat("Query1:", sqlString1, "\n")
    queryResult1 <- RMariaDB::dbGetQuery(impDbConn, sqlString1)
    thisConfidence <- 0
    if (nrow(queryResult1) == 0) {
      if (verbose) cat("FAILED - NO RESULT - TRY GENE2 FIELD\n")
      # first query failed, so try the other order, since it is an undirected
      # graph, and the database does not guarantee any ordering
      sqlString2 <- paste("select * from interactions where gene1='", gene2Name,
                          "' and gene2='", gene1Name, "'", sep = "")
      if (verbose) cat("Query2:", sqlString2, "\n")
      queryResult2 <- RMariaDB::dbGetQuery(impDbConn, sqlString2)
      if (nrow(queryResult2) == 0) {
        if (verbose) cat("FAILED - NO RESULT - SET TO ZERO\n")
        # second query failed, no hope; set to 0
        thisConfidence <- 0
        regainMatrix[gene1Idx, gene2Idx] <<- 0
        regainMatrix[gene2Idx, gene1Idx] <<- 0
      } else {
        if (verbose) cat("SUCCESS - PROCEED TO CONFIDENCE CHECK\n")
        # we know the second query was successful
        # now does it pass the confidence threshold?
        thisConfidence <- as.numeric(queryResult2[1, INTERACTIONS_TABLE_CONFIDENCE_FIELD_IDX])
        if (thisConfidence >= confidenceThreshold) {
          if (verbose) cat("Condfidence value:", thisConfidence,
                           "PASSES threshold:", confidenceThreshold, "\n")
          if (writeSifFile) {
            cat(gene1Name, gene2Name, thisConfidence, "\n", sep = "\t",
                file = sifFileCon, append = TRUE)
          }
        } else {
          if (verbose) cat("Confidence value:", thisConfidence, "DOES NOT PASS threshold:",
                           confidenceThreshold, ", setting matrix element to zero\n")
          thisConfidence <- 0
        }
        regainMatrix[gene1Idx, gene2Idx] <<- thisConfidence
        regainMatrix[gene2Idx, gene1Idx] <<- thisConfidence
      }
    } else {
      if (verbose) cat("SUCCESS\n")
      # we know the first query was successful
      thisConfidence <- as.numeric(queryResult1[1, INTERACTIONS_TABLE_CONFIDENCE_FIELD_IDX])
      # now does it pass the confidence threshold?
      if (thisConfidence >= confidenceThreshold) {
        if (verbose) cat("Confidence value:", thisConfidence, "PASSES threshold:",
                         confidenceThreshold, "\n")
        # yes it does!
        if (writeSifFile) {
          cat(gene1Name, gene2Name, thisConfidence, "\n",
              file = sifFileCon, sep = "\t", append = TRUE)
        }
      } else {
        if (verbose) cat("Confidence value:", thisConfidence, "DOES NOT PASS threshold:",
                         confidenceThreshold, "\n")
        thisConfidence <- 0
      }
      regainMatrix[gene1Idx, gene2Idx] <<- thisConfidence
      regainMatrix[gene2Idx, gene1Idx] <<- thisConfidence
    }
    data.frame(gene1Idx = gene1Name, gene2Idx = gene2Name, confidence = thisConfidence)
  })
  resultsDF <- do.call(rbind, resultsList)
  if (!is.null(sifFileCon) && writeSifFile) {
    close(sifFileCon)
  }

  list(filledMatrix = regainMatrix,
       filledGeneList = geneList,
       sif = resultsDF[resultsDF[, INTERACTIONS_TABLE_CONFIDENCE_FIELD_IDX] != 0, ])
}

computeProportionsFromCounts <- function(impDbConn = NULL,
                                         inGenesList = NULL,
                                         confidenceThreshold = 0,
                                         outGenesListTargetSize = 0,
                                         verbose = FALSE) {
  if (is.null(impDbConn)) {
    stop("impDbConn, a required parameter, is NULL")
  }
  if (is.null(inGenesList)) {
    stop("inGenesList, a required parameter, is NULL")
  }
  inGeneListSize <- length(inGenesList)
  if (outGenesListTargetSize < length(inGenesList)) {
    stop("outGenesListTargetSize", outGenesListTargetSize,
         "is less the length of the inGenesList", length(inGenesList))
  }

  # ap: select gene1, count(*) as count FROM interactions WHERE weight > 0.6 and
  # gene1 in('ZNF433','RHOXF2') group by gene1 order by count asc;
  # we need to assign a proportion to each gene in the initial list
  # a gene that has a larger proportion will get more genes added than one
  # with a smaller proportion
  if (verbose) cat("-------------------------------------------------------------\n")
  proportionalCounts <- vector(mode = "integer", length = inGeneListSize)
  for (geneFieldName in c('gene1', 'gene2')) {
    query = paste("select ", geneFieldName, ", count(*) as count FROM interactions ",
                  "WHERE weight >= ", confidenceThreshold, " AND ", geneFieldName, " ",
                  "IN( ", paste("'", inGenesList, "'", sep = "", collapse = ", "), " ) ",
                  "GROUP by ", geneFieldName, sep = "")
    if (verbose) cat(query, "\n")
    queryResult <- RMariaDB::dbGetQuery(impDbConn, query)
    # compute counts
    for (rowIdx in seq_along(rownames(queryResult))) {
      if (verbose) cat(rowIdx, inGenesList[rowIdx], "has", queryResult[rowIdx, 2], "in",
                       geneFieldName, "\n")
      proportionalCounts[rowIdx] <- proportionalCounts[rowIdx] + queryResult[rowIdx, 2]
    }
  }
  if (verbose) cat("Counts found: ", length(proportionalCounts), "\n")
  proportionArray = proportionalCounts / sum(proportionalCounts)
  # update true counts to use to form reGAIN
  proportionSum <- sum(proportionArray)
  proportionalCounts <- round(proportionArray * outGenesListTargetSize)
  # adjust the count array if doesn't total up to target size due to round remainders
  countSum <- sum(proportionalCounts)
  if (countSum != outGenesListTargetSize) {
    if (verbose) cat("countSum: ", countSum, " != outGenesListTargetSize: ", outGenesListTargetSize,
                     ", adjusting\n")
    remainder <- abs(outGenesListTargetSize - countSum)
    if (verbose) cat("Adjusting counts to remainder ", remainder, "positive or negative\n")
    subtractInLoop <- ifelse(remainder < 0, TRUE, FALSE)
    # use random indices?
    # randomCountIndices <- runif(n = remainder, min = 1, max = inGeneListSize)
    # NO, pick from the top genes and adjust the counts
    indicesDF <- data.frame(geneIdx = seq_along(inputGeneList), prop = proportionArray)
    rankedIndices <- indicesDF[order(indicesDF$prop, decreasing = TRUE), ]
    topIdx <- 1
    while (remainder && (topIdx < nrow(rankedIndices))) {
      if (verbose) cat("\t******************************************************\n")
      if (verbose) cat("Remainder ", remainder, "\n")
      if (verbose) cat(proportionalCounts, "\t")
      # pick a gene count index and adjust
      thisTopGeneIdx <- rankedIndices[topIdx, 1]
      if (subtractInLoop) {
        if (proportionalCounts[thisTopGeneIdx] > 0) {
          proportionalCounts[thisTopGeneIdx] <- proportionalCounts[thisTopGeneIdx] - 1
        }
      } else {
        proportionalCounts[thisTopGeneIdx] <- proportionalCounts[thisTopGeneIdx] + 1
      }
      # TODO: kill this section using random indices
      # randomCountIdx <- randomCountIndices[remainder]
      # if (subtractInLoop) {
      #   if (proportionalCounts[randomCountIdx] > 0) {
      #     proportionalCounts[randomCountIdx] <- proportionalCounts[randomCountIdx] - 1
      #   }
      # } else {
      #   proportionalCounts[randomCountIdx] <- proportionalCounts[randomCountIdx] + 1
      # }
      topIdx <- topIdx + 1
      remainder <- remainder - 1
    }
  }
  if (verbose) cat("proportionArray:", proportionArray, "*", outGenesListTargetSize, "=",
                   (proportionArray * outGenesListTargetSize), "\n")
  if (verbose) cat("Proportion sum:", proportionSum, "\n")
  if (verbose) cat("Counts to use: ", proportionalCounts, "\n")
  if (verbose) cat("Counts sum:", countSum, "\n")

  # sanity checks from constructing count arrays
  if (proportionSum < 0.999) {
    stop("Sum of proportions should equal 1 but its value is: ", proportionSum, "\n")
  }
  countSum <- sum(proportionalCounts)
  if (countSum < outGenesListTargetSize) {
    stop("Counts sum: ", countSum, " should equal: ", outGenesListTargetSize, "\n")
  }

  list(counts = proportionalCounts, probs = proportionArray)
}

addGenesWithCounts <- function(impDbConn = NULL,
                               inGenesList = NULL,
                               confidenceThreshold = 0,
                               proportionalCounts = NULL,
                               outGenesListTargetSize = 0,
                               verbose = FALSE) {
  if (is.null(impDbConn)) {
    stop("impDbConn, a required parameter, is NULL")
  }
  if (is.null(inGenesList)) {
    stop("inGenesList, a required parameter, is NULL")
  }
  inGeneListSize <- length(inGenesList)
  if (outGenesListTargetSize < inGeneListSize) {
    stop("outGenesListTargetSize", outGenesListTargetSize,
         "is less the length of the inGenesList", inGeneListSize)
  }
  if (is.null(proportionalCounts)) {
    stop("proportionalCounts, a required parameter, is NULL")
  }

  # ap: now for each gene in our genelist, add new nodes according to the counts
  # stored in proportionalCounts here I am querying both colums of the IMP database
  # with the query gene, ing the results to their respective counts
  # obtained proportionally, sorting the results and retaining the top "count"
  # of the gene; bcw: correct for cumulative rounding error to integers
  if (verbose) cat("-------------------------------------------------------------\n")
  dbGeneFieldNames <- c("gene1", "gene2")
  allAddedGenes <- inGenesList
  potentialNewGenes <- NULL
  keepAdding <- TRUE
  retries <- 0
  while (keepAdding) {
    for (countIdx in seq_along(proportionalCounts)) {
      thisGene <- inGenesList[countIdx]
      thisGeneAddCount <- proportionalCounts[countIdx]
      if (thisGeneAddCount > 0) {
        if (verbose) cat("-------------------------------------------------------------\n")
        if (verbose) cat(thisGene, "=>", thisGeneAddCount, "number of genes to add\n")
        # add the associated genes for each gene1 and gene2 fields
        potentialNewGenes <- NULL
        for (queryFieldIdx in seq_along(dbGeneFieldNames)) {
          geneFieldName <- dbGeneFieldNames[queryFieldIdx]
          otherFieldIdx <- ifelse(queryFieldIdx == 1, 2, 1)
          querySuccess <- FALSE
          genesWereAdded <- FALSE
          query <- paste("SELECT * FROM interactions WHERE " ,
                         geneFieldName , "='" , thisGene , "'" ,
                         " AND weight >= ", confidenceThreshold,
                         " ORDER BY weight DESC", sep = "")
          queryResult <- RMariaDB::dbGetQuery(impDbConn, query)
          queryResultLength <- nrow(queryResult)
          if (queryResultLength > 0) {
            # query success
            if (verbose) cat("Query success:\n", query, "\n")
            countDown <- thisGeneAddCount
            rowIdx <- 1
            while (countDown && (rowIdx <= length(queryResultLength))) {
              if (verbose) cat("Found:", queryResultLength,
                               "countDown:", countDown,
                               "rowIdx:", rowIdx, "\n")
              potentialGene <- queryResult[rowIdx, otherFieldIdx]
              potentialConfidence <- queryResult[rowIdx, INTERACTIONS_TABLE_CONFIDENCE_FIELD_IDX]
              if (verbose) cat("potentialGene:", potentialGene,
                               "potentialConfidence:", potentialConfidence, "\n")
              if ((sum(potentialGene == inGenesList) == 0) &&
                  (sum(potentialGene == allAddedGenes) == 0)) {
                potentialNewGenes <- rbind(potentialNewGenes,
                                           data.frame(gene = potentialGene,
                                                      confidence = potentialConfidence))
                countDown <- countDown - 1
                genesWereAdded <- TRUE
              } else {
                if (verbose) cat("Duplicate gene:", potentialGene, potentialConfidence, "skipped\n")
              }
              rowIdx <- rowIdx + 1
            }
            if (is.null(potentialNewGenes)) {
              querySuccess <- FALSE
            } else {
              querySuccess <- TRUE
              if (rowIdx > length(queryResultLength)) {
                if (verbose) cat("Exhausted query results, added", (thisGeneAddCount - countDown),
                                 "rows\n")
              }
              if (countDown == 0) {
                if (verbose) cat("Found the exact number of genes", thisGeneAddCount, "to add\n")
              }
            }
          } else {
            if (verbose) cat(thisGene, "Query did not find any results:\n", query, "\n")
          }

          if (genesWereAdded) {
            if (verbose) cat("Query SUCCESS, found", nrow(potentialNewGenes), "potentialNewGenes\n")
            finalGenesToAdd <- c()
            addStatus <- "Attempt to add new genes from potentialNewGenes table"
            if (verbose) cat("Sorting", nrow(potentialNewGenes),
                             "potentital genes by confidence to pick:",
                             thisGeneAddCount, "\n")
            allQueryGenesSorted <- as.character(potentialNewGenes[order(potentialNewGenes[, 2],
                                                                        decreasing = TRUE), 1])
            thisGeneAddCountAdj <- thisGeneAddCount
            if (thisGeneAddCount >= length(allQueryGenesSorted)) {
              thisGeneAddCountAdj <- length(allQueryGenesSorted)
              addStatus <- paste("Gene(s) added with proportional count adjusted from",
                                 thisGeneAddCount, "to QUERY RESULTS LENGTH",
                                 length(allQueryGenesSorted))
            }
            if ((length(allAddedGenes) + thisGeneAddCount) > outGenesListTargetSize) {
              thisGeneAddCountAdj <- outGenesListTargetSize - length(allAddedGenes)
              addStatus <- paste("Gene(s) added with proportional count adjusted from",
                                 thisGeneAddCount, "to the difference:",
                                 (length(allAddedGenes) - outGenesListTargetSize))
            }
            # add the genes found
            if (thisGeneAddCountAdj > 0) {
              sortedTopGenes = allQueryGenesSorted[1:min(thisGeneAddCountAdj,
                                                         length(allQueryGenesSorted))]
              finalGenesToAdd <- c(finalGenesToAdd, sortedTopGenes)
              addStatus <- paste("Gene(s) added with exact count, ",
                                 "or truncated to query size found",
                                 length(sortedTopGenes))
            } else {
              if (verbose) cat(thisGene, "thisGeneAddCountAdj < 1:",
                               thisGeneAddCountAdj, "\n")
            }
            if (verbose) cat("Gene add status:", addStatus, "\n")

            allAddedGenes <- unique(c(allAddedGenes, finalGenesToAdd))
            if (verbose) cat("New genes to add:", length(finalGenesToAdd),
                             finalGenesToAdd,"\n")
            if (verbose) cat("Genes added so far:", length(allAddedGenes), "\n")
          } else {
            warning("Query successful but no genes were added for: ", thisGene)
          }
        }
      }
    }

    # TODO: probably don't need this block any more???
    if (length(unique(allAddedGenes)) == outGenesListTargetSize) {
      cat("Number of genes added reached the target, WE OUT\n")
      keepAdding <- FALSE
    } else {
      cat("Attempting retry\n")
      retries <- retries + 1
      if (retries > RETRY_LIMIT) {
        cat("Number of retries exceeds the :", RETRY_LIMIT, "BAILING\n")
        keepAdding <- FALSE
      } else {
        # retry using one more gene in the  part of the SQL statement
        cat("Increasing the count array and running retry:", retries, "RETRY\n")
        cat("********************** R E T R Y ******************************\n")
        proportionalCounts <- proportionalCounts + retries
        keepAdding <- TRUE
      }
    }
  }

  allAddedGenes
}

# ap: this function handles adding interacting genes to the network from the IMP
# database i.e. add more user specified number of connections meeting a certain
# threshold - bcw - 9/6/17 - check for duplicates to insure desirec counts
addConnectionsFromDatabase <- function(impDbConn = NULL,
                                       geneList = NULL,
                                       scoreList = NULL,
                                       hasMainEffectScores = FALSE,
                                       ranked = FALSE,
                                       confidenceThreshold = 0,
                                       targetListSize = 0,
                                       sifFilename = "",
                                       verbose = FALSE) {
  if (verbose) cat("-------------------------------------------------------------\n")
  # sane parameters?
  if (is.null(impDbConn)) {
    stop("impDbConn is a required parameter")
  }
  if (is.null(geneList)) {
    stop("geneList is a required parameter")
  }
  inGeneListSize <- length(geneList)
  cat("addConnectionsFromDatabase: the original gene list has", inGeneListSize, " genes\n")
  if (targetListSize < 1) {
    stop("targetListSize must be > 0, or use fillRegain()")
  }

  if (verbose) cat("-------------------------------------------------------------\n")
  cat("Computing gene proportions from queried input gene counts\n")
  proportionResults <- computeProportionsFromCounts(impDbConn = impDbConn,
                                                    inGenesList = geneList,
                                                    confidenceThreshold = confidenceThreshold,
                                                    outGenesListTargetSize = targetListSize,
                                                    verbose = verbose)
  print(proportionResults$counts)
  if (verbose) cat("-------------------------------------------------------------\n")
  cat("Adding genes based on probabilistically computed counts\n")
  addGenesUniqueResults <- addGenesWithCounts(impDbConn = impDbConn,
                                              inGenesList = geneList,
                                              confidenceThreshold = confidenceThreshold,
                                              proportionalCounts = proportionResults$counts,
                                              outGenesListTargetSize = targetListSize,
                                              verbose = verbose)
  print(addGenesUniqueResults)

  if (verbose) cat("-------------------------------------------------------------\n")
  if (verbose) cat("Original genes:", length(geneList), ":", geneList, "\n\n")
  if (verbose) cat("Added genes:", length(addGenesUniqueResults),
                   ":", addGenesUniqueResults, "\n\n")
  newGeneList <- unique(c(geneList, addGenesUniqueResults))
  if (verbose) cat("New gene list, unique(original + added), size:", length(newGeneList),
                   ":", newGeneList, "\n\n")

  if (verbose) cat("-------------------------------------------------------------\n")
  if (verbose) cat("Running fillRegain()\n")
  fillRegain(impDbConn = impDbConn,
             geneList = newGeneList,
             scoreList = scoreList,
             hasMainEffectScores = hasMainEffectScores,
             ranked = ranked,
             confidenceThreshold = confidenceThreshold,
             addMode = FALSE,
             sifFilename = "",
             verbose = FALSE)
}
