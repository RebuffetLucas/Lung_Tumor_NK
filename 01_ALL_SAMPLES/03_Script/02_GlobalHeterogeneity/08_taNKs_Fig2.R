library(Seurat)
library(EnhancedVolcano)
library(ComplexHeatmap)
library(dplyr)
library(RColorBrewer)
library(escape) #v1.12.0
library(tibble)      
library(tidyr)    
library(circlize)
library(enrichplot)
library(clusterProfiler)
library(dittoSeq)
set.seed(1234)

NK <- readRDS("path/to/NK_only.RDS")

# -- DimPlot - taNKs, non taNKs -- Fig. S3E

NK$tumor <- ifelse(NK$Fine_Annotation %in% c("NK3_CAMK4", "NK3_ENTPD1"), 
                    "tumor-associated", "non tumor-associated")

NK$tumor[NK$Fine_Annotation == "ILC1"] <- "ILC1"
NK$tumor[NK$Fine_Annotation == "ILC3"] <- "ILC3"
NK$tumor <- factor(NK$tumor, levels = c("tumor-associated",  "non tumor-associated", 
                                        "ILC1","ILC3"))


Colors <- c("#6699FF","#E15759" ,"#999999", "#CCCCCC")

pdf(file = "~/Multiome/analysis/Objects_plots/plots/FINAL/NK_Tumor/Final/NK_FINAL_subsetted/NormalizeData/v5/WNN/tumor_unique/Dimplot_tur.pdf",
    width = 3, height = 3)
DimPlot(NK, reduction = "WNN_UMAP", label = T, group.by = "tumor", cols = Colors,
        label.size = 4,label.box = T, repel = T, raster = F) & NoAxes() & NoLegend()
dev.off()

# -- Volcano - DEGs - taNKs vs. non taNKs -- Fig. 2F
NK <- SetIdent(NK, value = NK$tumor)

DEGs_taNKs <- FindMarkers(NK, assay = "RNA", ident.1 = "tumor-associated", ident.2 = "non tumor-associated",logfc.threshold = 0,
                          min.pct = 0.25, test.use = "MAST")

pdf("path/to/Volcano_taNK.pdf", height = 6, width = 7)
EnhancedVolcano(DEGs_taNKs,
                lab = DEGs_taNKs$X,
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = "tumor-resident lung NK cells",
                titleLabSize = 25,
                subtitleLabSize = 18,
                FCcutoff = 0.5,
                ylab = bquote(~-Log[10] ~ italic(P[adj])),
                xlim = c(-4, 4),
                pointSize = 3.0,
                legendIconSize = 3.0,
                labSize = 4,
                pCutoff = 0.0001,
                legendPosition = "bottom",
                selectLab = c("ENTPD1", "GZMA", "ITGAE", "ITGA1", "RGS1", "LAYN", "CTLA4",
                              "CD96", "CXCR6", "IL32", "CCL5", "KLRC1", "KLRC2",
                              "PLPP1", "CSF1", "IRF4", "CCL4", "ITGA2", "CXCL13","PDCD1",
                              "LAG3","TOX","DAPK2","LINGO2",
                              "NKG7", "NFKBIA","FCGR3A", "CCL3", "CD69", "GZMM","PRF1",
                              "ZEB1", "CX3CR1","FGFBP2", "SPN2", "TNFRSF9"),
                drawConnectors = T,
                colConnectors = 'black',
                gridlines.minor = F, 
                gridlines.major = F,
                #labFace = 'bold',
                boxedLabels = TRUE, 
                colAlpha = 4/5,
                max.overlaps = Inf,
                col = c("#CCCCCC","#CCCCCC","#CCCCCC","#CC3333")) +
                theme(plot.title = element_blank(),
                      plot.subtitle = element_blank(),
                      axis.text = element_text(color = 'black'),
                      axis.title = element_text(color = 'black'),
                      axis.line = element_line(color = 'black'),
                      legend.text = element_blank(),
                      legend.position = "none")
dev.off()

# -- Signature Enrichments -- Fig. 2G

NK <- readRDS("path/to/NK_only.RDS")

#T cell signatures to test
dyf_Li <- readLines("path/to/dysfunctional_Li.txt") # Li et al. 2019 (doi: 10.1016/j.cell.2018.11.043.)
exh_JA <- readLines("path/to/CD8_exhausted_Jerby-Arnon.txt") # Jerby-Arnon et al. 2018 (https://doi.org/10.1016/j.cell.2018.09.006)
cyt_JA <- readLines("path/to/CD8_cytotoxic_Jerby-Arnon.txt") # Jerby-Arnon et al. 2018 (https://doi.org/10.1016/j.cell.2018.09.006)
termexh_Oli <- readLines("path/to/term_ex_Oli.txt") # Oliveira et al. 2021 (https://doi.org/10.1038/s41586-021-03704-y)
tum_reac <- readLines("path/to/Lowery_tumor-reactive_CD8.txt") # Lowery et al. 2022 (DOI: 10.1126/science.abl5447)
tum_reac_top30 <- tum_reac[1:30]

trNK <- readLines("path/to/trNK_Marquardt.txt") # Marquardt et al. 2019 (doi: 10.1038/s41467-019-11632-9)
adapt <- readRDS("path/to/Signature_adaptive_Rückert.rds")
top30_adaptiveNK <- rownames(adapt[order(-adapt$avg_log2FC), ])[1:30]

#prepare the Zheng et al. 2021 (DOI: 10.1126/science.abe6474)

CD8_Zheng <- read.csv("path/to/Pan_cancer_Zheng_CD8.csv")
CD8_Zheng %>%
  group_by(cluster.name) %>%
  slice_max(n = 30, order_by = comb.ES) -> top30_CD8_Zheng

CD8_Zheng_top30 = subset(top30_CD8_Zheng, select = c("geneSymbol","cluster.name"))

sig_CD8_Zheng <- split(CD8_Zheng_top30, CD8_Zheng_top30$cluster.name)
Trm_Zheng <- sig_CD8_Zheng$`CD8.c10(ZNF683+CXCR6+ Trm)`
Ttex_Zheng <- sig_CD8_Zheng$`CD8.c12(terminal Tex)`
Trm_Zheng <- as.vector(Trm_Zheng$geneSymbol)
Ttex_Zheng <- as.vector(Ttex_Zheng$geneSymbol)

sig <- list(dysfunctional_Li=dyf_Li, exhauted_Jerby_Arnon=exh_JA, cytotoxic_Jerby_Arnon=cyt_JA, 
            term_exh_Oliveira=termexh_Oli, tumor_reactive_Lowery =tum_reac_top30,
            Trm_Zheng=Trm_Zheng, Ttex_Zheng=Ttex_Zheng, trNK=trNK, top30_adaptiveNK=top30_adaptiveNK)

#enrich signature in taNKs vs. non taNKs
agg_taNK <- AggregateExpression(NK, assays = "RNA",return.seurat = T,
                                group.by = c("orig.ident", "tumor"))

agg_taNK_noseu <- AggregateExpression(NK, assays = "RNA",return.seurat = F,
                                       group.by = c("orig.ident", "tumor"))


agg_taNK@meta.data
agg_taNK@meta.data[c("Patient", "cluster")] <- do.call(rbind, strsplit(agg_taNK$orig.ident,"_"))

agg_taNK$cluster <- factor(agg_taNK$cluster, levels = c("tumor-associated", "non tumor-associated", 
                                                          "ILC1", "ILC3"))

ES <- enrichIt(obj = as.matrix(agg_taNK_noseu[["RNA"]]),
               gene.sets = sig, 
               groups = 1000, cores = 12, 
               min.size = 5, method = "UCell")

agg_taNK <- AddMetaData(agg_taNK, ES)
ES_taNK <- data.frame(agg_taNK[[]])

#without ILCs
ES_taNK <- subset(ES_taNK, tumor == "tumor-associated"| tumor == "non tumor-associated")

#Heatmap
#get a matrix and scale
mat_plotting <- ES_taNK %>%
  dplyr::select(-orig.ident, -tumor, -Patient, -cluster) %>%
  scale() %>%
  t() %>%
  as.matrix()

mypalette <- rev(brewer.pal(n = 9, "RdBu"))
myBreaks <- seq(-2, 2, length.out = 9)

#set colnames manually
colnames(mat_plotting)<- c(
  "dysfunctional signature",
  "exhausted signature",
  "cytotoxic signature",
  "tumor-reactive CD8 signature",
  "resident-memory CD8 signature",
  "terminal exhaustion CD8 signature",
  "tissue-resident NK signature",
  "adaptive NK signature")

ES_taNK$cluster <- factor(ES_taNK$cluster, levels = c("non tumor-associated", "tumor-associated"))
ES_taNK$cluster <- droplevels(ES_taNK$cluster) #ILCs were still saved as levels

h <- Heatmap(mat_plotting, name = "UCell \n z-score",  
             column_split = ES_taNK$cluster,
             cluster_columns = F,
             show_column_dend = FALSE,
             cluster_column_slices = T,
             column_title_gp = gpar(fontsize = 12),
             column_gap = unit(1, "mm"),
             cluster_rows =F,
             show_row_dend = FALSE,
             col=colorRamp2(myBreaks, mypalette),
             row_names_gp = gpar(fontsize = 12),
             column_title_rot = 0,
             show_column_names = F,
             use_raster = F,
             raster_quality = 4,
             row_names_side = "left")

pdf("path/to/HeatMap_taNK.pdf", width = 8, height = 3)
draw(h, padding = unit(c(2, 4, 2, 2), "mm")) 
dev.off()


# -- Signature Enrichments -GSEA -- Fig. S3F

#order DEGs of taNKs vs. non taNKs by log2FC
DEGs_taNKs <- DEGs_taNKs[order(DEGs_taNKs$avg_log2FC, decreasing = T),]

#make a named vector
gl <-DEGs_taNKs$avg_log2FC
names(gl)<-make.names(DEGs_taNKs$X, unique = T)

#read in signature to test
tum_reac <- read.delim("path/to/Lowery_tumor-reactive_CD8.txt", header = FALSE) %>%
  rename(gene = V1) %>%
  mutate(term = "CD8_tumor-reactive_Lowery") %>%
  select(term, gene)

#perform enrichment
GSEA_tumreac <- GSEA(gl, TERM2GENE = tum_reac,
                     pvalueCutoff = 0.05, seed=T)

#saveResults
write.csv(GSEA_tumreac, "path/to/GSEA_result_tumor_reactive_Lowery.csv")

pdf("path/to/GSEA_tumor-reactive.pdf",height = 4, width = 7)
gseaplot2(GSEA_tumreac, geneSetID = "CD8_tumor-reactive_Lowery", base_size = 10.5)
dev.off()

# -- Volcano - DEGs - NK3_ENTPD1 vs. NK3_CAMK4 -- Fig. S3G

NK <- SetIdent(NK, value = NK$Fine_Annotation)
DEGs_CAMKvsENTPD1 <- FindMarkers(NK, assay = "RNA", ident.1 = "NK_ENTPD1", ident.2 = "NK_CAMK4",logfc.threshold = 0,
                                   test.use = "MAST")

write.csv(DEGs_CAMKvsENTPD1, "path/to/ENTPD1vsCAMK4_DEGs.csv")

pdf("path/to/Volcano_ENTvsCAMK.pdf",height = 6, width = 7)
EnhancedVolcano(DEGs_CAMKvsENTPD1,
                lab = DEGs_CAMKvsENTPD1$X,
                x = 'avg_log2FC',
                y = 'p_val_adj',
                titleLabSize = 25,
                subtitleLabSize = 18,
                FCcutoff = 0.5,
                ylab = bquote(~-Log[10] ~ italic(P[adj])),
                pointSize = 3.0,
                legendIconSize = 3.0,
                xlim = c(-4, 4),
                labSize = 4,
                pCutoff = 0.0001,
                legendPosition = "bottom",
                selectLab = c("ZBTB16", "GZMA", "PLPP1", "CRTAM", "KLRC1", "KLRB1",
                              "GZMH", "FCER1G", "TMIGD2", "ENTPD1", "CD44", "KLRF2",
                              "ITGAE", "IRF4", "EOMES", "CXCL13", "A2M", "CD247", "CCL5",
                              "KIR2DL4", "ZNF831", "BCL11B", "NFATC2", "TNFRSF9", "INPP4B"),
                drawConnectors = T,
                colConnectors = 'black',
                gridlines.minor = F, 
                gridlines.major = F,
                boxedLabels = TRUE, 
                colAlpha = 4/5,
                max.overlaps = Inf,
                col = c("#CCCCCC","#CCCCCC","#CCCCCC","#CC3333")) +
                theme(plot.title = element_blank(),
                      plot.subtitle = element_blank(),
                      axis.text = element_text(color = 'black'),
                      axis.title = element_text(color = 'black'),
                      axis.line = element_line(color = 'black'),
                      legend.text = element_blank(),
                      legend.position = "none")
dev.off()

# -- Signature Enrichment - adaptive NKs in NK3_ENTPD1 vs. NK3_CAMK4 -- Fig. S3H

adapt <- readRDS("path/to/Signature_adaptive_Rückert.rds")
top30_adaptiveNK <- rownames(adapt[order(-adapt$avg_log2FC), ])[1:30]


#aggregate
agg <- AggregateExpression(subset(NK, Fine_Annotation == "NK3_CAMK4"|Fine_Annotation == "NK3_ENTPD1"), 
                            assays = "RNA",return.seurat = T,
                            group.by = c("orig.ident", "Fine_Annotation"))

agg_noseu <- AggregateExpression(subset(NK, Fine_Annotation == "NK3_CAMK4"|Fine_Annotation == "NK3_ENTPD1"), 
                                 assays = "RNA",return.seurat = F,
                                 group.by = c("orig.ident", "Fine_Annotation"))


agg@meta.data
agg$Fine_Annotation <- gsub('-', '_', agg$Fine_Annotation)

ES <- enrichIt(obj = as.matrix(agg_noseu[["RNA"]]),
               gene.sets = list(top30_adaptive=top30_adaptiveNK), 
               groups = 1000, cores = 12, 
               min.size = 5, method = "UCell")

agg <- AddMetaData(agg, ES)
ES_taNK <- data.frame(agg[[]])

#BoxPlot using DittoSeq
counts <- ES %>% select(top30_adaptive) %>% as.matrix() %>% t()
metadata <- ES %>% select(-top30_adaptive) %>% as.matrix()

bulkSCE <- importDittoBulk(
  x = list(counts = counts),
  metadata = data.frame(metadata))

Colors <- c("#FF994E","#6699FF")

pdf("path/to/Boxplot_adaptiveNK.pdf", 
dittoBoxPlot(bulkSCE, "top30_adaptive", group.by = "Fine_Annotation", color.panel = Colors, 
             boxplot.width = 0.5,
             boxplot.color = "black",boxplot.lineweight = 0.3, min = 0,
             ylab = "adaptive NK signature score", jitter.size = 1.2) &
  geom_pwc(group.by="x.var",label="p.signif",hide.ns=TRUE, method = "wilcox_test",
           vjust = 0.3) &
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
        axis.text.y = element_text(size = 8,color = "black"),
        axis.text.x = element_text(size = 12,  color = "black", angle = 45),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, color = "black"),
        legend.position = "none")
dev.off()










