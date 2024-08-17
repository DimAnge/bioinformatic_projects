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
library(tidyverse)
library(GEOquery)
library(biomaRt)


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

#Annotation using biomart
mart <- useEnsembl(biomart = 'ensembl', 
                   dataset = 'hsapiens_gene_ensembl', 
                   version = 110)

#Query based on ensemble_gene_id
# Remove the ".number" suffix from row names so i can query
annot_res <- res_clean
rownames(annot_res) <- sub("\\.\\d+$", "", rownames(annot_res))
annot <- getBM(filters= "ensembl_gene_id", #which identifier you are using /ex ENSG00053523
               attributes= c("ensembl_gene_id", #what you retrieve
                             "description", 
                             "start_position", 
                             "end_position", 
                             "strand", 
                             "hgnc_symbol"), #which attributes you want to collect
               values= rownames(annot_res), #the names of your genes
               mart= mart)
annot_res$ensembl_gene_id <- rownames(annot_res)
res_df <- as.data.frame(annot_res)
final <-  res_df %>%
  left_join(annot, by = 'ensembl_gene_id') 

#Exporting annotated results
write.xlsx(final, file="data/DE_Annotated_Results_GSE260628.xlsx", colNames=TRUE)

#Creating plots
#MA plot 
MA_plot <-plotMA(res)

#Adding significance level (logFC 0.58 (x1.5 differentially expressed) padj < 0.5)
final_plot <- final
final_plot <- final_plot %>%
  mutate(significance = case_when(
    log2FoldChange >= 0.58 & padj <= 0.05 ~ 'Up-regulated',
    log2FoldChange <= -0.58 & padj <= 0.05 ~ 'Down-regulated',
    abs(log2FoldChange) < 0.58 | padj > 0.05 ~ 'Not significant'))

#Volcano plot 
volcano <- ggplot(data=final_plot, 
                  aes(x=log2FoldChange, y=-log10(padj),
                      color = significance,)) + 
  geom_point() + 
  geom_vline(xintercept=c(-0.58, 0.58), col="blue") +  #logFC threshold
  geom_hline(yintercept=-log10(0.05), col="blue") + #padj threshold
  labs(x = 'LogFC', y = '-Log10(PValue-adjusted)', color = 'Significance') +
  scale_color_manual(values = c('blue', 'grey', 'red')) +
  theme_minimal()

#Exporting plot
pdf("Plots/Volcano_plot_GSE260628.pdf", height = 8, width = 8)
print(volcano)
dev.off() 

#Heatmap
#I will select the top 25 Up and top 25 Downn-regulated genes
up25 <- final_plot %>%
  slice_max(order_by = log2FoldChange, n = 25)
down25 <- final_plot %>%
  slice_min(order_by = log2FoldChange, n = 25)
top50 <- bind_rows(up25,down25)

#For the heatmap we will use the rlog transformation
rlog <- rlog(dds, blind = FALSE)
rownames(rlog) <- sub("\\.\\d+$", "", rownames(rlog)) #covert rownames to format needed
top50_rlog <- assay(rlog)[top50$ensembl_gene_id, ]
#Convert ensembl_gene_id names to hgnc symbols 
gene_to_symbol <- setNames(annot$hgnc_symbol, annot$ensembl_gene_id)
new_rownames <- gene_to_symbol[rownames(top50_rlog)]
rownames(top50_rlog) <- new_rownames

heatmap <- Heatmap(top50_rlog,
                   row_labels = rownames(top50_rlog), 
                   row_names_gp = gpar(fontsize = 5), 
                   column_names_gp = gpar(fontsize = 7), 
                   heatmap_legend_param = list(title = "Rlog\nexpression"), 
                   top_annotation = HeatmapAnnotation(Condition = meta$treatment, 
                                                      which = 'column', 
                                                      col = list(Condition = c(
                                                        "untreated" = 'green',
                                                        "heparin lyase" = 'red')
                                                      )        )        )

#Exporting Heatmap
pdf("Plots/Heatmap_GSE260628.pdf", height = 8, width = 8)
print(heatmap)
dev.off()

#Session info
sink(file = 'data/session_info_GSE260628_DE_analysis_DESeq2.txt')
sessionInfo()
sink()

