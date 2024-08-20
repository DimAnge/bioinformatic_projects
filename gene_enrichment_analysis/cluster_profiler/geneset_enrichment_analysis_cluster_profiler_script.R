#Dimitris Angelakis - Aug 2024
#Gene set enrichment analysis of RNA-seq data of human lung microvascular ECs using ClusterProfiler
#In this mini project i will use public data to perform a DEA using the R packages : ClusterProfiler for practice.
#Overrepresentation analysis (ORA), Functional Class Scoring (FCS)
#Data : Richter RP, Odum JD, Margaroli C, Cardenas JC et al. Trauma promotes heparan sulfate modifications and cleavage that disrupt homeostatic gene expression in microvascular endothelial cells. Front Cell Dev Biol 2024;12:1390794. PMID: 39114570
#Data available at Gene Expression Omnibus (GEO) database (Acession Number GSE260628), on NCBI

#Load required libraries
library(openxlsx)
library(ggplot2)
library(dplyr)
library(org.Hs.eg.db)
library(DOSE)
library(enrichplot)
library(clusterProfiler)
library(msigdbr)
library(GSEABase)
library(tibble)

#Reading DE results
degs <- read.xlsx('data/DE_Annotated_Results_GSE260628.xlsx')

#Thresholds
logFC_threshold = 0.58 
padj_threshold=0.05   
#Annotate
degs <- degs %>%
  mutate(regulation = case_when(
    log2FoldChange > logFC_threshold & padj < padj_threshold ~ 'Upregulated',
    log2FoldChange < -logFC_threshold & padj < padj_threshold ~ 'Downregulated',
    .default = 'Not significant'
  ))

#Checks volcano plot
ggplot(degs, aes(x = log2FoldChange,
                 y=-log10(padj),          
                 color = regulation)) + 
  geom_point() +
  theme_minimal() + 
  geom_hline(yintercept = -log10(padj_threshold)) + 
  geom_vline(xintercept = logFC_threshold) + 
  geom_vline(xintercept = -logFC_threshold)

#Make variable for the availiable up and down regulated genes
upgenes <- degs %>%
  filter(regulation == 'Upregulated', 
         !is.na(hgnc_symbol),   
         !duplicated(hgnc_symbol)) %>%
  pull(hgnc_symbol)             

downgenes <- degs %>%
  filter(regulation == 'Downregulated', 
         !is.na(hgnc_symbol),
         !duplicated(hgnc_symbol)) %>%
  pull(hgnc_symbol)

#All genes in the DEG dataset that have available hgnc symbols
background <- degs %>%
  filter(!is.na(hgnc_symbol)) %>%
  filter(!duplicated(hgnc_symbol)) %>%
  pull(hgnc_symbol)

#Gathering pathway gene sets 
#Molecular Signature Database

genesets1 <- msigdbr(species = 'Homo sapiens',
                     category= 'H') %>%         #a   Hallmark pathways
  bind_rows(msigdbr(species = 'Homo sapiens',
                    category= 'C2',
                    subcategory = 'CP'))%>%   #α canonical pathways
  bind_rows(msigdbr(species = 'Homo sapiens',
                    category= 'C5',
                    subcategory = 'GO:BP'))   #a gene ontology:biological process

#getting term sizes
#gsea package default -> filtered >10 genes
term_sizes <- genesets1 %>%
  distinct(gs_name,gene_symbol) %>%  #unique kai ta dio
  group_by(gs_name) %>%           #group by pathway
  summarise(count = n()) 

#Filter genesets (>10)
sets_to_keep <- term_sizes %>%
  filter(count > 10) %>%
  pull(gs_name)

genesets <- genesets1 %>%
  filter(gs_name %in% sets_to_keep)

#clusterProfiler - Over Expression Analysis
#Term to genes to filter out any non unique combinations
t2g_cp <- genesets %>%
  dplyr::distinct(gs_name,gene_symbol) %>%  
  as.data.frame()
#Up-regulated genes
ora_cp_up <- enricher(upgenes,
                      pvalueCutoff = 0.05,
                      pAdjustMethod = 'BH', #same as fdr/padj
                      universe = background,
                      qvalueCutoff = 0.2,
                      TERM2GENE = t2g_cp) 
res_cp_ora_up <- ora_cp_up@result %>%
  filter(p.adjust < 0.05) 
#Down-regulated genes
ora_cp_down <- enricher(downgenes,
                        pvalueCutoff = 0.05,
                        pAdjustMethod = 'BH',
                        universe = background,
                        qvalueCutoff = 0.2,
                        TERM2GENE = t2g_cp)

res_cp_ora_down <- ora_cp_down@result %>%
  filter(p.adjust < 0.05)

#Plotting clusterProfiler ORA
cPr_plot <- res_cp_ora_up %>%
  mutate(regulation = 'Upregulation') %>% 
  bind_rows(res_cp_ora_down %>%
              mutate(regulation = 'Downregulation')) %>%
  mutate(source = gsub('^([[:alpha:]]+)_[[:print:]]+','\\1',ID)) %>% 
  mutate(ID = gsub('_',' ',gsub('^[[:alpha:]]+_','',ID)))            
#Converting GeneRatio to numeric
temp <- as.numeric(unlist(strsplit(cPr_plot$GeneRatio, split = '/')))

cPr_plot$GeneRatio2 <- temp[c(TRUE, FALSE)]/temp[c(FALSE, TRUE)]

top_cPr_plot <- cPr_plot %>%
  slice_max(order_by = GeneRatio2, n = 20)       

#Plotting
pdf('Plots/clustPr_ora_bubble_plot.pdf', width = 7, height = 7)
ggplot(top_cPr_plot,
       aes(x = GeneRatio2,
           y = reorder(ID,GeneRatio2),
           color = regulation,
           size = -log10(p.adjust))) + 
  geom_point() + 
  scale_y_discrete(breaks = cPr_plot$ID,
                   labels = strtrim(cPr_plot$ID,
                                    width = 30)) + 
  labs(x = 'Gene ratio',
       y = 'Term',
       color = 'Up/Downregulation',
       size = '-Log10(Adj. PValue)') + 
  scale_color_manual(breaks = c('Downregulation','Upregulation'),
                     values = c('lightblue',
                                'red4')) + 
  facet_wrap(~ source,
             scales = 'free_y')
dev.off()

#FCS -Functional Class Scoring (FCS)- ranked gene list
#Creates a ranked named vector based on DE results 
ranked <- degs %>%
  mutate(rank_metric = stat) %>%   
  arrange(desc(rank_metric)) %>%                    
  filter(!is.na(hgnc_symbol), hgnc_symbol != '') %>%
  select(hgnc_symbol, rank_metric)

#Change to the format needed in gsea of clusterProfiler
ranked_vector <- ranked$rank_metric            
names(ranked_vector) <- ranked$hgnc_symbol

#FCS using GSEA of clusterProfiler
gsea_cp <- GSEA(geneList = ranked_vector,
                pvalueCutoff = 0.05,
                pAdjustMethod = 'BH',
                eps = 0,             
                seed = 149,          #reproducibility seed
                minGSSize = 1,      #limits for geneset sizes
                maxGSSize = 1000,
                TERM2GENE = t2g_cp)  #provide genesets 

#gsea result object so we use '@'
res_fcs_cp <- gsea_cp@result

#Plotting clusterProfiler FCS 
#We use the enrichment score to see if a path is up or down regulated
cP_fcs_plot <- res_fcs_cp %>%
  mutate(source = gsub('^([[:alpha:]]+)_[[:print:]]+','\\1',ID),
         regulation = ifelse(NES > 0, 'Upregulation','Downregulation')) %>%    #if paths normalized enrichment score>0 -> upregulated otherwise downregulated
  mutate(ID = gsub('_',' ',gsub('^[[:alpha:]]+_','',ID))) %>%
  filter(source == 'GOBP') %>%
  slice_max(order_by = abs(NES), n = 40)                    

# Plotting bubble for FCS
pdf('Plots/clustPr_FCS_bubble_plot.pdf', width = 13, height = 7)
ggplot(cP_fcs_plot,
       aes(x = NES,
           y = reorder(ID, NES),
           color = regulation,
           size = -log10(p.adjust))) + 
  geom_point() + 
  scale_y_discrete(breaks = cP_fcs_plot$ID,
                   labels = strtrim(cP_fcs_plot$ID,
                                    width = 50)) + 
  labs(x = 'NES',
       y = 'Term',
       color = 'Up/Downregulation',
       size = '-Log10(Adj. PValue)') + 
  scale_color_manual(breaks = c('Downregulation', 'Upregulation'),
                     values = c('lightblue', 'red4')) 
dev.off()

#Barplot
pdf('Plots/clustPr_FCS_barplot_plot.pdf', width = 15, height = 10)
ggplot(cP_fcs_plot,
       aes(x = NES,
           y = reorder(ID,NES),
           fill = regulation,
           label = p.adjust)) + 
  geom_bar(stat = 'identity') + 
  scale_y_discrete(breaks = cPr_plot$ID,
                   labels = strtrim(cPr_plot$ID,
                                    width = 40)) + 
  labs(x = 'NES',
       y = 'Term',
       fill = 'Up/Downregulation') + 
  scale_fill_manual(breaks = c('Downregulation','Upregulation'),
                    values = c('lightblue',
                               'red4')) + 
  theme_minimal() + 
  facet_wrap(~ source,
             scales = 'free_y')
dev.off()

#Session info
sink('data/session_info_enrichment_analysis_cp.txt')
sessionInfo()
sink()
