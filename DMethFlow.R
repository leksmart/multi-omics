# Set root and main directories
root <- "/mmfs1/home/david.adeleke/Thesis"
main <- "DNAMeth"
proj.scope <- c("PAAD-US", "PACA-CA", "PACA-AU")
# Create directory if it doesn't exist and set working directory
uscontrol= c("DO32751", "DO218973", "DO50272", "DO32769")

if (!dir.exists((file.path(root, main)))) {dir.create((file.path(root, main,  recursive = TRUE)))}

setwd(file.path(root, main))

# Load necessary libraries
required_packages <- c("InfiniumPurify","circlize", "tidyverse", "ComplexHeatmap", 
                       "SummarizedExperiment", "TCGAbiolinks", "ggplot2", "DT", "tidyr",
                       "reshape2", "plyr", "plotly", "tibble", "hrbrthemes", "viridis",
                       "GGally", "kableExtra", "survival", "survminer", "grid", "gridExtra",
                       "rmdformats","R.utils","readr","data.table", "edgeR", "limma",
                       "org.Hs.eg.db", "DESeq2","enrichplot", "DOSE", "forcats", "ggstance",
                       "ReactomePA", "clusterProfiler", "GOSemSim", "magrittr", "minfi", "IlluminaHumanMethylation450kanno.ilmn12.hg19",
                       "RColorBrewer", "missMethyl", "minfiData", "Gviz", "DMRcate","stringr")




# Load libraries, install if not available
# Load libraries, install if not available
lapply(required_packages, function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
})


coho=1

source("script.R")




