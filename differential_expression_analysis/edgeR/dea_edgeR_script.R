#Dimitris Angelakis - Aug 2024
#Differential expression analysis of RNA-seq data of human bladder cancer cells using edgeR
#In this mini project i will use public data to perform a DEA using the R package edgeR for practise.
#Data : Wenjie W (2024), Human bladder cancer cell line high throughput sequencing - genes regulated by hMSH2 - untreated. Gene Expression Omnibus (GEO) database (Accession Number GSE193754), Available on NCBI.

##Loading libraries
library(edgeR)
library(ComplexHeatmap)
library(circlize)
library(openxlsx)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(FactoMineR)
library(GEOquery)

##Reads expression
expr1 <- read.table(file = "data/GSE193754_All_reads_counts.txt", header = TRUE, row.names = 1)
expr<-expr1[,c("EJ_Cas9_1","EJ_Cas9_2","EJ_Cas9_3","EJ_KO_1","EJ_KO_2","EJ_KO_3")]

##Reading series matrix (sample information)
gse<-getGEO(filename = "data/GSE193754_series_matrix.txt",GSEMatrix = TRUE, getGPL = FALSE)
metadata <- pData(gse)
meta <- metadata %>%
  select(c("title", "genotype:ch1")) %>%
  rename(sample_title = title, sample_type = "genotype:ch1") 

##Data cleaning, preparation and normalization
#Rearranges data and checks if order is identical between expr data and metadata
rownames(meta) <- NULL
meta <- meta %>%
  arrange(match(sample_title,c("EJ_Cas9_1", "EJ_Cas9_2", "EJ_Cas9_3", "EJ_KO_1", "EJ_KO_2", "EJ_KO_3"))) %>%
  column_to_rownames(var="sample_title")
identical(colnames(expr),rownames(meta))

#Checks if expr data contain NAs
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


##Performing a Principal Component Analysis with FactorMineR (just for practice/small sample number)
#Transposing the table of logcpm
trans_log_cpm <- t(log_cpm)

#PCA (keeping the scale.unit since we have removed most lowly expressed genes)
pca_cpm <- PCA(trans_log_cpm,scale.unit = T, graph = F)

#Extract the coordinates from the PCA data object
pca_cord <- as.data.frame(pca_cpm$ind$coord)

#Variance explained % by each component in pca_cpm$eig
pca_perc <- pca_cpm$eig[,2]

#Plotting the PCA with ggplot2 using the first 2 Principal Components
pca_plot <- ggplot(pca_cord, aes(pca_cord[,1], pca_cord[,2],
       color = meta$sample_type))+
  geom_point(key_glyph = "point") +
  labs(title = "PCA - GSE193754",
       x = paste0("PC1(", round(pca_perc[1], 2), "%)"),
       y = paste0("PC2(", round(pca_perc[2], 2), "%)"),
       color = 'Sample Type') +
  theme(plot.title = element_text(hjust = 0.5))
print(pca_plot)

#Exporting plot
pdf(file = 'Plots/PCA_GSE193754.pdf', width = 8, height = 6)
print(pca_plot)
dev.off()

##Differential Expression Analysis
#Estimating Dispersion (in edgeR : GLMCommonDisp, GLMTrendedDisp, GLMTagwiseDisp )
dd1 <- estimateDisp(dd, design, verbose = TRUE)

#Fitting the model (estimates the effect of each variable and evaluate the logFC)
fit <- glmQLFit(dd1, design)

#Comparison between  MSH2-deficient vs wild type  (create contrast)
colnames(design) <- c("grouphMSH2_deficient","groupwild_type")
contr <- makeContrasts(grouphMSH2_deficient - groupwild_type, levels = design)

#Running DE test
deg_lrt <- glmQLFTest(fit, contrast = contr)
head(deg_lrt$table)

##False Discovery Rate Correction (Benjamin - Hochberg) - Multiple comparison adjustment
#Adding fdr
deg_lrt$table$fdr <- p.adjust(deg_lrt$table$PValue, method = "BH")

##Filtering results
#fdr < 0.05 and absolute lofFC >= 0.58 (strict 0.01 and 1)
top <- deg_lrt$table[c(deg_lrt$table$fdr <= 0.05) & abs(deg_lrt$table$logFC) >= 0.58,] #34 genes

#Exporting Results
write.xlsx(top, file = "data/DE_Results_GSE193754_FDR_0_05_logFC_0_58.xlsx", rowNames = TRUE, colNames = TRUE)

#Saving workspace
save.image("data/GSE193754_DE_results.Rdata")

##Visualizing results
#Adding significance level at the data frame
res <- deg_lrt$table
res <- res %>%
  mutate(significance = case_when(
    logFC >= 0.58 & fdr <= 0.05 ~ 'Up-regulated',
    logFC <= -0.58 & fdr <= 0.05 ~ 'Down-regulated',
    abs(logFC) < 0.58 | fdr > 0.05 ~ 'Not significant'))

#Volcano plot (top right significant up-regulated, top left significant down-regulated)
  volcano <- ggplot(data=res, 
         aes(x=logFC, y=-log10(fdr),
         color = significance,)) + 
    geom_point() + 
    geom_vline(xintercept=c(-0.58, 0.58), col="blue") +  #logFC threshold
    geom_hline(yintercept=-log10(0.05), col="blue") + #fdr threshold
    labs(x = 'LogFC', y = '-Log10(FDR)', color = 'Significance') +
    scale_color_manual(values = c('blue', 'grey', 'red')) +
    theme_minimal()
    
#Exporting plot
pdf("Plots/Volcano_plot_GSE193754.pdf", height = 8, width = 8)
print(volcano)
dev.off()    

#Heatmap  (We use the logCPM value) of the 34 differential expressed genes
top_log_cpm <- log_cpm[rownames(top),]

#Creating the heatmap
heatmap <- Heatmap(top_log_cpm,
        row_labels = rownames(top_log_cpm), 
        row_names_gp = gpar(fontsize = 5), 
        column_names_gp = gpar(fontsize = 7), 
        heatmap_legend_param = list(title = "LogCPM\nexpression"), 
        top_annotation = HeatmapAnnotation(Condition = meta$sample_type, 
                                           which = 'column', 
                                           col = list(Condition = c(
                                                                    "wild type" = 'green',
                                                                    "hMSH2-deficient" = 'red')
                                           )        )        )

#Exporting Heatmap
pdf("Plots/Heatmap_GSE193754.pdf", height = 8, width = 8)
print(heatmap)
dev.off()

#Session info
sink(file = 'data/session_info_GSE193754_DE_analysis_edgeR.txt')
sessionInfo()
sink()
