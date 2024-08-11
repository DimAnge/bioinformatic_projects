#Dimitris Angelakis - Aug 2024
#Differential expression analysis of RNA-seq data of human bladder cancer cells using edgeR
#In this mini project i will use public data to perform a DEA using the R package edgeR for practise.
#Data : Wenjie W (2024), Human bladder cancer cell line high throughput sequencing - genes regulated by hMSH2 - untreated. Gene Expression Omnibus (GEO) database (Accession Number GSE193754), Available on NCBI.

##Loading libraries
library(edgeR)
library(biomaRt)
library(ComplexHeatmap)
library(circlize)
library(openxlsx)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(FactoMineR)
library(plotly)
library(GEOquery)

##reads expression
expr1 <- read.table(file = "data/GSE193754_All_reads_counts.txt", header = TRUE, row.names = 1)
expr<-expr1[,c("EJ_Cas9_1","EJ_Cas9_2","EJ_Cas9_3","EJ_KO_1","EJ_KO_2","EJ_KO_3")]

##reading matrix (sample information)
gse<-getGEO(filename = "data/GSE193754_series_matrix.txt",GSEMatrix = TRUE, getGPL = FALSE)
metadata <- pData(gse)
meta <- metadata %>%
  select(c("title", "genotype:ch1")) %>%
  rename(sample_title = title, sample_type = "genotype:ch1") 

##data cleaning, preparation and normalization
#rearranges data and checks if order is identical betweenn expr data and metadata
rownames(meta) <- NULL
meta <- meta %>%
  arrange(match(sample_title,c("EJ_Cas9_1", "EJ_Cas9_2", "EJ_Cas9_3", "EJ_KO_1", "EJ_KO_2", "EJ_KO_3"))) %>%
  column_to_rownames(var="sample_title")
identical(colnames(expr),rownames(metadata_sub))

#checks if expr data contain NAs
table(is.na(expr))

#Determine factor (wild type vs hMSH2-deficient)
factor(meta$sample_type)
group <- as.factor(meta$sample_type)

#Creating a model matrix
design <- model.matrix(~group+0)
design

#Creating a differential gene expression object with edgeR
d <- DGEList(counts=expr, group = group)

#Calculating the normalization factors (Need to be used before the cpm() function, which returns normalized counts)
d <- calcNormFactors(d)

##Filtering lowly expressed genes
#keeping genes with > 1 cpm for at least 5% of the samples (since we have only 6 sample,we just keep genes with cpm > 1)
keep <- rowSums(cpm(d) > 1) >= ceiling(0.05*dim(d)[2])
d_new <- d[keep,]
dim(d_new)[1] # 16786 genes were kept from the 56303
(dim(d_new)[1]/dim(d)[1])*100 #or ~=29.81% genes were kept

##Creating a new DGEList object
dd <- DGEList(counts = d_new, group = group)
#Recalculating factors
dd <- calcNormFactors(dd)
#Normalizing counts with cpm()
norm_d <- cpm(dd)
#Getting the log values of  our data
log_cpm <- cpm(dd, log=TRUE)

#Exporting normalized data
write.table(norm_d, file ="data/GSE193754_cpm.txt", sep = '\t')
write.table(log_cpm, file ="data/GSE193754_logcpm.txt", sep = '\t')


##Performing a Principal Component Analysis

