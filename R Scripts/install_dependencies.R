## script checks for and install dependencies

check.packages <- function(pkg){
  # check.packages function: install and load multiple R packages.
  # Check to see if packages are installed. Install them if they are not, 
  # then load them into the R session.
  # https://gist.github.com/smithdanielle/9913897
  
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, repos = "http://cran.us.r-project.org", dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

if (!("devtools" %in% installed.packages()[,"Package"])){
  install.packages("devtools", repos = "http://cran.us.r-project.org", dependencies = TRUE)
}
library(devtools)

# install github package
if (!("Rinbix" %in% installed.packages()[,"Package"])){
  devtools::install_github("insilico/Rinbix", build_vignettes = TRUE)
}

# install cran packages
packages <- c("ggplot2", "CORElearn", "gridExtra", "lattice", "caret", "glmnet", "xgboost", "igraph",
              "gdata", "SDMtools")

# install bioconductor packages
source("https://bioconductor.org/biocLite.R")
biocLite("GEOquery")
biocLite("Biobase")
biocLite("preprocessCore")
