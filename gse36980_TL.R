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
library(DESeq2)
library(devtools)
library(goseq)
library(oligo)
#PCA
library(simpleaffy)
library(tidyverse)
library(ggfortify)
#oligo
library(RCurl)
library(biomaRt)

#Metadata import

celFiles <- list.files("GSE36980_RAW_TL/", full =TRUE)
gse <- read.celfiles(celFiles)
gse36980 <- getGEO(filename = "GSE36980_series_matrix.txt")
metadata <- gse36980@phenoData@data
metadata <- metadata[34:62,]
CN <- metadata[c("title")]
CN_1 <- metadata[c("title")]

#Data normalisation

rma <- rma(gse) 

#PCA plot

colour <- c(rep('AD-TL', 10), rep('Normal-TL', 19))
rawGSE <- exprs(rma)
pca_OG <- prcomp(t(rawGSE), scale. = T, center = T)
pca_OG_df <- as.data.frame(pca_OG$x)
PCA_raw <- ggplot(pca_OG_df, aes(x=PC1, y=PC2, color=colour))+ geom_point() + stat_ellipse() +
  labs(title = "Raw Data PCA Plot")


#rma normalisation boxplot 
df <- exprs(rma)
rma_boxplot <- boxplot(df, main="Relative Signal BoxPlot map", ylab="Relative log expression signal-RMA", las=2)

#Annotation
ID <- rownames(df)
mart <- useEnsembl(biomart = "ensembl", 
                   dataset = "hsapiens_gene_ensembl", 
                   mirror = "useast")
gene_list <- getBM(attributes=c('affy_hugene_1_0_st_v1', 'hgnc_symbol'), 
                   filters = 'affy_hugene_1_0_st_v1', 
                   values = ID, 
                   mart = mart)
gene_list_df <- as.data.frame(gene_list)
duplicate_probeID <- gene_list_df[!duplicated(gene_list_df$hgnc_symbol),]
e <- exprs(rma)
table_merge <- merge(x = duplicate_probeID, y = e, by.x = "affy_hugene_1_0_st_v1", by.y = "row.names")
table_merge <- table_merge[!duplicated(table_merge$hgnc_symbol),]
table_merge <- na.omit(table_merge)
rownames(table_merge) <- table_merge$hgnc_symbol
annotated <- table_merge[-c(1,2)]

#Gene filtering

mean <- rowMeans(annotated)
remove_lower_0.02 <- annotated[which(mean > 0.02),]
geneFilt <- quantile(mean, p=0.02)
geneFilt <- annotated[which(mean>geneFilt),]

#limma
#creating model + volcano plot
ADorNor <- data.frame(patient = metadata$title)
rownames(ADorNor) <- rownames(metadata)
ADorNor$patient <- ifelse(str_detect(ADorNor$patient, regex(".AD")), 1, 0)
matrix <- model.matrix(~ patient, ADorNor)
fit <- limma::lmFit(remove_lower_0.02, matrix)
efit <- eBayes(fit)
lod <- -log10(efit[["p.value"]][,2])
mtstat<- efit[["t"]][,2]
limma_output <- topTable(efit, coef = 2, adjust = "fdr", n = 50000)
volcano <- EnhancedVolcano(toptable = limma_output, lab = rownames(limma_output), x = "logFC", y = "P.Value", pCutoff = 0.05, FCcutoff = 0.5)

#heatmap
group <- data.frame(patient = metadata$title)
top50 <- topTable(efit, number=50) 
filt_50 <- data.frame(geneFilt[rownames(geneFilt) %in% rownames(top50),])
heatmap_4107 <- pheatmap(filt_50, main="Top 50 DEG", 
                         labels_row = rownames(filt_50),  
                         labels_col = colnames(filt_50),
                         annotation_row=group, 
                         annotation_col=group, 
                         annotation_colors=list(patient=c(AD ="red", Normal="blue")))

#DEG 

top10 <- rownames(topTable(efit, n = 10))

#MODULE 5
#DEG 
logFC <- limma_output %>% dplyr::select('logFC')
logFC_vec <- as.vector(logFC$logFC)
element_names <- rownames(limma_output)
names(logFC_vec) <- element_names

#threshold + filtering (sorted, named, numeric vector)
filtered_FC <- logFC_vec[logFC_vec > 0.4]
arrange_FC <- sort(filtered_FC, decreasing = TRUE)
topDEG <- data.frame(arrange_FC)

#selecting symbol and entrezid
Entrezid_symbol <- getBM(
  attributes = c('hgnc_symbol', 'entrezgene_id'),
  filters = 'affy_hugene_1_0_st_v1',
  values = ID,
  mart = mart
)
df_entrezid <- Entrezid_symbol %>% dplyr::select('hgnc_symbol', "entrezgene_id")
df_entrezid <- na.omit(df_entrezid)
names <- names(arrange_FC)
selected <- df_entrezid[df_entrezid$hgnc_symbol %in% names, ]
selected <- unique(selected)
rownames(selected) <- 1:nrow(selected)

#gene ontology + boxplots
CC <- enrichGO(selected$entrezgene_id, OrgDb = org.Hs.eg.db, ont = "CC" , readable = T )
barplot_GO_CC <- barplot(CC)

MF <- enrichGO(selected$entrezgene_id, OrgDb = org.Hs.eg.db, ont = "MF" , readable = T )
barplot_GO_MF <- barplot(MF)

BP <- enrichGO(selected$entrezgene_id, OrgDb = org.Hs.eg.db, ont = "BP" , readable = T )
barplot_GO_BP <- barplot(BP)

#KEGG
KEGG <- enrichKEGG(selected$entrezgene_id, pvalueCutoff = 0.05)
dotplot_KEGG <- dotplot(KEGG)

#Gene-concept network
convert <- setReadable(KEGG, OrgDb=org.Hs.eg.db, keyType = "ENTREZID")
cnetplot_GCN <- cnetplot(convert, foldChange = arrange_FC, categorySize = "pvalue")

#Global/universal gene set enrichment analysis (GSEA)
gene_set <- msigdbr(species = "Homo sapiens", category = "H")
h <- gene_set %>% select (gs_name, entrez_gene)

removed_duplicate <- Entrezid_symbol[!duplicated(Entrezid_symbol$SYMBOL),]
removed_duplicate <- na.omit(removed_duplicate)
rowname_symbol <- removed_duplicate$SYMBOL
rownames(removed_duplicate) <- rowname_symbol
GSEA_merge <- merge(x = removed_duplicate, y = limma_output, by = "row.names")
GSEA_logFC <- as.vector(GSEA_merge$logFC)
GSEA_entrezid <- as.vector(GSEA_merge$ENTREZID)
names(GSEA_logFC) <- GSEA_entrezid
sorted <- sort(GSEA_logFC, decreasing = TRUE)

GSEA_analysis <- GSEA(sorted, TERM2GENE = h)
GSEA_plot <- gseaplot(GSEA_analysis, geneSetID = 2)

for_module_6 <- names(filtered_FC)
write.table(arrange_FC, file = 'DEG.txt')

