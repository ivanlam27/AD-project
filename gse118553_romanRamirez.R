# Roman Ramirez, rr8rk@virginia.edu
# Transcriptomics Analysis on GSE118553

library(WGCNA)
library(EnhancedVolcano)
library(ggplot2)
library(affyQCReport)
library(GEOquery)
library(simpleaffy)
library(tidyverse)
library(hgu133plus2.db)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(limma)
library(pheatmap)
library(edgeR)
library(Biobase)
library(affyQCReport)
library(affy)
library(simpleaffy)
library(affyPLM)
library(GEOquery)
library(utilsIPEA)
library(clusterProfiler)
library(enrichplot)
library(magrittr)
library(tidyr)
library(ggnewscale)
library(msigdbr)
library(readxl)

#PCA
library(simpleaffy)
library(tidyverse)
library(ggfortify)

#metadata

gse <- ReadAffy(celfile.path = "data_gse118553/GSE118553_RAW")
series_matrix <- getGEO(filename = "data_gse118553/GSE118553_series_matrix.txt")
metadata <- series_matrix@phenoData@data
# CN <- metadata[C("title")]