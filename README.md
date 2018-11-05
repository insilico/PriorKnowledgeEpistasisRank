# PriorKnowledgeEpistasisRank
Epistasis network centralities that incorporate prior knowledge
Author: Saeid Parvandeh and Brett McKinney

### How to install required packages
    > check.packages <- function(pkg){
      # check.packages function: install and load multiple R packages.
      # Check to see if packages are installed. Install them if they are not, 
      # then load them into the R session.
      # https://gist.github.com/smithdanielle/9913897
  
      new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
      if (length(new.pkg)) 
       install.packages(new.pkg, repos = "http://cran.us.r-project.org", dependencies = TRUE)
      sapply(pkg, require, character.only = TRUE)
    > }

    > if (!("devtools" %in% installed.packages()[,"Package"])){
        install.packages("devtools", repos = "http://cran.us.r-project.org", dependencies = TRUE)
    > }
    > library(devtools)

    # install github package
    > if (!("Rinbix" %in% installed.packages()[,"Package"])){
        devtools::install_github("insilico/Rinbix", build_vignettes = TRUE)
    > }

    # install cran packages
    > packages <- c("ggplot2", "CORElearn", "gridExtra", "lattice", "caret", "glmnet", "xgboost", "igraph",
              "gdata", "SDMtools")

    # install bioconductor packages
    > source("https://bioconductor.org/biocLite.R")
    > biocLite("GEOquery")
    > biocLite("Biobase")
    > biocLite("preprocessCore")
    

### How to find vignettes
You will be able to find two vignettes, pathway.Rmd and classification.Rmd. The first one manifests the steps of creating pathaway enrichment and how prior knowledge improve Reactome Pathways for MDD, and the second shows the classification accuracies where we add prior knowledge. 

### Accuracy plot
Training accuracy (Cambridge data) and independent validation accuracy (Japan data) with centrality feature selection without prior knowledge (left panels) and with prior knowledge (right panels). Top: co-expression net-work centrality feature selection methods, PageRank (PR) and Katz. Bottom row: expression-epistasis network centrality methods, EpistasisRank (ER) and EpistasisKatz (EK). Accuracies computed by xgboosted trees with nested cross-validation. Xgboost accuracies without feature selection also shown (squares).

![Accuracy plots](Acc_original_plot.png)