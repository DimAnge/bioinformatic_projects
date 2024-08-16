#Dimitris Angelakis - Aug 2024
#Differential expression analysis of RNA-seq data of human lung microvascular ECs using DESeq2
#In this mini project i will use public data to perform a DEA using the R package DESeq2 for practice.
#Data : Richter RP, Odum JD, Margaroli C, Cardenas JC et al. Trauma promotes heparan sulfate modifications and cleavage that disrupt homeostatic gene expression in microvascular endothelial cells. Front Cell Dev Biol 2024;12:1390794. PMID: 39114570
#Data available at Gene Expression Omnibus (GEO) database (Acession Number GSE260628), on NCBI


##Loading libraries
library(DESeq2)
library(ComplexHeatmap)
library(circlize)
library(openxlsx)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(FactoMineR)
library(GEOquery)


##Reads expression
expr <- read.xlsx("data/GSE260628_counts.xlsx")
expr <- column_to_rownames(expr, var = "EnsemblID")

##Reading series matrix (sample information)
gse<-getGEO(filename = "data/GSE260628_series_matrix.txt",GSEMatrix = TRUE, getGPL = FALSE)
metadata <- pData(gse)
meta <- metadata %>%
  select(c("title", "treatment:ch1")) %>%
  rename(sample_title = title,  treatment = "treatment:ch1")
rownames(meta) <- NULL

##Data cleaning, preparation and normalization
#Rearranges data and checks if order is identical between expr data and meta
meta <- meta %>%
  arrange(match(sample_title,c(colnames(expr)))) %>%
  column_to_rownames(var="sample_title")
identical(colnames(expr),rownames(meta))

#Creating a DESeq dataset object
dds <- DESeqDataSetFromMatrix(countData = expr,
                       colData = meta,
                       design = ~ treatment)

#Filtering lowly expressed genes (<10 reads) from 67065 genes to 19620
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#Set factor level, untreated as reference
dds$treatment <- relevel(dds$treatment, ref = "untreated")

#Run DESeq
dds <- DESeq(dds)
res <- results(dds)

#Results
res_0_01 <- results(dds, alpha = 0.01)
summary(res_0_01)

##Filtering top results
#padj < 0.05 and absolute lofFC >= 0.58 (strict 0.01 and 1)
res_clean <- na.omit(res)
top <- res_clean[c(res_clean$padj <= 0.05) & abs(res_clean$log2FoldChange) >= 0.58,] #679 genes
top$EnsemblID <- rownames(top)
top <- top[, c("EnsemblID", setdiff(names(top), "EnsemblID"))]

#Exporting Results
write.xlsx(top, file = "data/DE_Results_GSE260628_Padj_0_05_logFC_0_58.xlsx", rowNames = FALSE, colNames = TRUE)

#Creating plots
#MA plot (with adj Pvalue < 0.05)
plotMA(res)
