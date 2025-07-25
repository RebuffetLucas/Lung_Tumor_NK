library(Seurat)
library(Signac)
library(ggsignif)
library(ggpubr)
library(dittoSeq)
library(escape) #v1.12.0
library(tidyverse)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(rlist)
library(reshape2)
set.seed(1234)

NK <- readRDS("path/to/NK_only.RDS")

# -- Gene Activity analysis --

DefaultAssay(NK) <- "peaks"

gene.activities <- GeneActivity(NK)

NK[['GeneActivity']] <- CreateAssayObject(counts = gene.activities)

NK <- NormalizeData(
  object = NK,
  assay = 'GeneActivity',
  normalization.method = 'LogNormalize',
  scale.factor = median(NK$nCount_ATAC)
)

# -- Marker calculations Coarse and Fine --

#calculate Marker genes for coarse clusters (res 0.2) and fine clusters (res 0.4)

#Coarse 0.2

NK <- SetIdent(NK, value = NK$wsnn_res.0.2)
Markers_coarse <- FindAllMarkers(NK, assay = "RNA", min.pct = 0.25, test.use = "MAST")

write.csv(Markers_coarse, "path/to/Marker_coarse.csv")

#top30 per cluster
Markers_coarse  %>%
  group_by(cluster) %>%
  top_n(n = 30, wt = avg_log2FC) -> top30_coarse

write.csv(top30_coarse, "path/to/Marker_coarse_top30.csv")

#Fine 0.4 

NK <- SetIdent(NK, value = NK$wsnn_res.0.4)
Markers_fine <- FindAllMarkers(NK, assay = "RNA", min.pct = 0.25, test.use = "MAST")

write.csv(Markers_fine, "path/to/Marker_fine.csv")

#top30 per cluster
Markers_fine  %>%
  group_by(cluster) %>%
  top_n(n = 30, wt = avg_log2FC) -> top30_fine

write.csv(top30_fine, "path/to/Marker_fine_top30.csv")

# -- DotPlots -- using genes from the top20 markers/known NK cell genes that describe subsets 

#Coarse 0.2 - RNA - Fig. S2A

pdf(file = "path/to/DotPlot_coarse_RNA.pdf", width =5.5, height = 2)
DotPlot(NK, group.by = "wsnn_res.0.2",
        features =  c("NKG7", "FCGR3A", "FGFBP2", "S1PR5","CX3CR1","NCAM1", 
                      "SELL","XCL1", "CD44", "IL32", "CD3E", "CD52","ENTPD1", "ITGAE", "ITGA1",
                      "IL7R","RORA", "CD5", "CD6", "KIT", "AHR", "RORC"),
        assay= "RNA", cols = c("RdBu")) + RotatedAxis() +
  theme(axis.text.y = element_text(size = 9),
        axis.text.x = element_text(size = 9, face = "bold"),
        legend.title = element_text(size=9),
        legend.text = element_text(size=8),
        legend.key.size = unit(0.3, "cm"),
        axis.title = element_blank())
dev.off()

#Coarse 0.2 - Gene Activity - Fig. S2A

pdf(file = "path/to/DotPlot_coarse_GA.pdf", width =5, height = 2)
DotPlot(NK, group.by = "wsnn_res.0.2",
        features =  c("CXCR2", "B3GAT1", "CX3CR1", "ZNF462", "ZNF667", "CTLA4", "ITGA2", "CXCL13", "CCR6",
                      "IL7R", "CD3G", "B3GALT5", "RORC", "KIT"),
        assay= "GeneActivity", cols = c("PiYG")) + RotatedAxis() +
  theme(axis.text.y = element_text(size = 9),
        axis.text.x = element_text(size = 9, face = "bold"),
        legend.title = element_text(size=9),
        legend.text = element_text(size=8),
        legend.key.size = unit(0.3, "cm"),
        axis.title = element_blank())
dev.off()

#Fine 0.4 - RNA - Fig. 1F

pdf(file = "path/to/DotPlot_Fine_RNA.pdf", width =5.5, height = 2)
DotPlot(NK, group.by = "wsnn_res.0.4",
        features =  c("FGFBP2", "FCGR3A","CCL4", "CCL3","CX3CR1", "NFKB1",
                      "NCAM1","IL7R", "XCL1", "CD44","GZMK","CXCR4",
                      "IL32","CD3E", "CD52","ITGA1","CAMK4" ,"ENTPD1", "ITGAE", "GZMA",
                      "MKI67", "TOP2A","RORA", "CD5", "CD6", "KIT", "AHR", "RORC"),
        assay= "RNA", cols = c("RdBu")) + RotatedAxis() +
  theme(axis.text.y = element_text(size = 9),
        axis.text.x = element_text(size = 9, face = "bold"),
        legend.title = element_text(size=9),
        legend.text = element_text(size=8),
        legend.key.size = unit(0.3, "cm"),
        axis.title = element_blank())
dev.off()

#Fine 0.4 - Gene Activity - Fig. 1G

pdf(file = "path/to/DotPlot_Fine_GA.pdf", width =5, height = 2)
DotPlot(NK, group.by = "wsnn_res.0.4",
        features =  c("CXCR2", "B3GAT1", "CX3CR1", "ZNF462", "ZNF667", "CTLA4", "ITGA2", "CXCL13", "CCR6",
                      "IL7R", "CD3G", "B3GALT5", "RORC", "KIT"),
        assay= "GeneActivity", cols = c("PiYG")) + RotatedAxis() +
  theme(axis.text.y = element_text(size = 9),
        axis.text.x = element_text(size = 9, face = "bold"),
        legend.title = element_text(size=9),
        legend.text = element_text(size=8),
        legend.key.size = unit(0.3, "cm"),
        axis.title = element_blank())
dev.off()

# -- Annotation Based on Genes and Signature Scorings --
#Coarse Annotation/depends on clusters
NK$wsnn_res.0.2 <- as.numeric(NK$wsnn_res.0.2)

NK$Coarse_Annotation <- "NK3" 
NK$Coarse_Annotation[NK$wsnn_res.0.2 == 2] <- "NK1"
NK$Coarse_Annotation[NK$wsnn_res.0.2 == 3] <- "NK2"
NK$Coarse_Annotation[NK$wsnn_res.0.2 == 4] <- "ILC1"
NK$Coarse_Annotation[NK$wsnn_res.0.2 == 5] <- "ILC3"

NK$Coarse_Annotation <- factor(NK$Coarse_Annotation, 
                               levels = c("NK1", "NK2", "NK3", "ILC1", "ILC3"))

#Fine Annotation
NK$wsnn_res.0.4 <- as.numeric(NK$wsnn_res.0.4)

NK$Fine_Annotation <- "NK1_FGFBP2" 
NK$Fine_Annotation[NK$wsnn_res.0.4 == 2] <- "NK1_CCL4"
NK$Fine_Annotation[NK$wsnn_res.0.4 == 3] <- "NK_NFKB1"
NK$Fine_Annotation[NK$wsnn_res.0.4 == 4] <- "NK2"
NK$Fine_Annotation[NK$wsnn_res.0.4 == 5] <- "NK3_GZMK"
NK$Fine_Annotation[NK$wsnn_res.0.4 == 6] <- "NK3_CAMK4"
NK$Fine_Annotation[NK$wsnn_res.0.4 == 7] <- "NK3_ENTPD1"
NK$Fine_Annotation[NK$wsnn_res.0.4 == 8] <- "NKprol"
NK$Fine_Annotation[NK$wsnn_res.0.4 == 9] <- "ILC1"
NK$Fine_Annotation[NK$wsnn_res.0.4 == 10] <- "ILC3"

NK$Fine_Annotation <- factor(NK$Fine_Annotation, 
                             levels = c("NK1_FGFBP2", "NK1_CCL4", "NK_NFKB1", "NK2", 
                                        "NK3_GZMK", "NK3_CAMK4", "NK3_ENTPD1", "NKprol",
                                        "ILC1", "ILC3"))

# -- Verification using Signature Scorings -- BoxPlots + Heatmaps --

signatures <- readRDS("path/to/signatures_List.RDS") #List of signatures to test obtained from indicated papers
#13-gene NK signature, NK1, NK2, NK3, ILCvsNK, ILC1, ILC2, ILC3

#Coarse - aggregate Matrix + escape package v1.12.0 - (same done for fine clusters)

agg_Coarse <- AggregateExpression(NK, assays = "RNA",return.seurat = T, 
                                  slot = "scale.data", group.by = c("orig.ident", "Coarse_Annotation"))

agg_Coarse_noseu <- AggregateExpression(NK, assays = "RNA",return.seurat = F, slot = "scale.data",
                                        group.by = c("orig.ident", "Coarse_Annotation"))

agg_Coarse@meta.data[c("Patient", "cluster")] <- do.call(rbind, strsplit(agg_Coarse$orig.ident,"_"))

agg_Coarse$cluster <- factor(agg_Coarse$cluster, levels = levels(NK$Coarse_Annotation))

#enrich
ES_coarse <- enrichIt(obj = as.matrix(agg_Coarse_noseu[["RNA"]]), 
                      gene.sets = signatures,  
                      groups = 1000, cores = 12, 
                      min.size = 5, method = "UCell")

agg_Coarse <- AddMetaData(agg_Coarse, ES_coarse)
ES_df_coarse <- data.frame(agg_Coarse[[]])

#HeatMap - Fig. S2B
Colors_coarse <- c("#E15759","#86BCB6","#59A14F","#999999", "#CCCCCC") #color factors need to match annotation factor levels
Main_names <- names(signatures)[1:8] # select only the ones you want to plot (NK1-3, ILC sigs)
func <- names(signatures)[9:10] # phenotype (trNK, adaptive NK)

mat_NK_main <- ES_df_coarse %>% select(Main_names) %>% as.matrix()
mat_NK_main <- t(scale(mat_NK_main))

mypalette <- rev(brewer.pal(n = 9, "RdBu"))
myBreaks <- seq(-2, 2, length.out = 9)

pdf("path/to/Signatures_Sigs_Coarse.pdf", width = 7, height = 4)
Heatmap(mat_NK_main, name = "UCell \n z-score",  
        column_split = ES_df_coarse$cluster,
        cluster_columns = F,
        show_column_dend = FALSE,
        cluster_column_slices = T,
        column_title_gp = gpar(fontsize = 12),
        column_gap = unit(1, "mm"),
        cluster_rows =F,
        show_row_dend = FALSE,
        col=colorRamp2(myBreaks, mypalette),
        row_names_gp = gpar(fontsize = 12),
        column_title_rot = 45,
        top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = Colors_coarse))),
        show_column_names = F,
        use_raster = F,
        raster_quality = 4,
        row_names_side = "left")
dev.off()

#the same plotting was applied for Fine clusters, see Fig. S1C

#BoxPlot - Fig. 1C and D

#fake a dittobulk object for easy plotting
counts<- ES_df_coarse %>% select(Main_names, trNK_Marquardt, adaptive_Rückert) %>% as.matrix()
counts <- t(counts)
metadata <- ES_df_coarse %>% select(Patient, cluster)

bulkSCE_coarse <- importDittoBulk(
  x = list(counts = counts),
  metadata = data.frame(metadata))

plt1 <- list()
plt2 <- list()

max1 <- c(1.1, 1.4, 1.2,0.8, 1.1, 0.65, 0.2,1.1) #adjusted to y axis
max2 <- c(0.8, 0.65) #adjusted to y axis

y1 <- c("13-gene NK signature score","NK1 signature score", "NK2 signature score", "NK3 signature score", 
        "ILC signature score", "ILC1 signature score", "ILC2 signature score","ILC3 signature score")
y2 <- c("trNK signature score", "adaptive NK signature score")

step_increase1 <- c(0.18, 0.18, 0.18, 0.18, 0.18, 0.21, 0,0.21) # for plotting
step_increase2 <- c(0.25, 0.2)

#Main Signatures - Fid. 1C
for (i in 1:length(Main_names)){
  plt1[[i]] <- dittoBoxPlot(bulkSCE_coarse, Main_names[[i]], group.by = "cluster", color.panel = Colors_coarse, 
                            boxplot.width = 0.5,
                            boxplot.color = "black",boxplot.lineweight = 0.3, min = 0,
                            ylab = y1[[i]], main = y1[[i]], jitter.size = 1.2, max = max1[[i]]) &
    geom_pwc(group.by="x.var",label="p.adj.signif",hide.ns=TRUE, method = "dunn_test", step.increase = step_increase1[[i]]) & #multiples comp Dunns
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
          axis.text.y = element_text(size = 8,color = "black"),
          axis.text.x = element_text(size = 12,  color = "black", angle = 45),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 12, color = "black"),
          legend.position = "none")}

pdf("path/to/Boxplot_Coarse_NKsign.pdf", 
    width = 12, height = 5)
ggarrange(plotlist = plt1, nrow = 2,ncol = 4)
dev.off()

#NK phenotype signatures - Fid. 1D

for (i in 1:length(func)){
  plt2[[i]] <- dittoBoxPlot(bulkSCE_coarse, func[[i]], group.by = "cluster", color.panel = Colors_coarse, 
                            boxplot.width = 0.5,
                            boxplot.color = "black",boxplot.lineweight = 0.3, min = 0,
                            ylab = y2[[i]], main = y2[[i]], jitter.size = 1.2, max = max2[[i]]) &
    geom_pwc(group.by="x.var",label="p.adj.signif",hide.ns=TRUE, method = "dunn_test", step.increase = 0.1) &
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
          axis.text.y = element_text(size = 8,color = "black"),
          axis.text.x = element_text(size = 12,  color = "black", angle = 45),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 12, color = "black"),
          legend.position = "none")}

pdf("path/to/Boxplot_Coarse/Boxplot_Coarse_funcNKphenotype.pdf", 
    width =3, height = 5)
ggarrange(plotlist = plt2, nrow = 2,ncol = 1)
dev.off()

# --- plot histology and tumor stage per fine cluster

table <- dittoBarPlot(NK, var = "Fine_Annotation", group.by = "Donor_histo", retain.factor.levels = F,data.out=T)
table_stage <- dittoBarPlot(NK, var = "Fine_Annotation", group.by = "Donor_stage", retain.factor.levels = F,data.out=T)

cell_freq <- table$data

cell_freq_stage <- table_stage$data

#histo column
cell_freq$histo <- ifelse(cell_freq$grouping %in% c("BS840_squamous cell carcinoma", "BS1314_squamous cell carcinoma",
                                                    "BS1319_squamous cell carcinoma", "BS1322_squamous cell carcinoma"),
                          "squamous cell carcinoma", 
                          "adenocarcinoma")

cell_freq$histo <- factor(cell_freq$histo, levels=c("adenocarcinoma","squamous cell carcinoma"))

#stage column
cell_freq_stage$stage <- ifelse(cell_freq_stage$grouping %in% c("BS824_II", "BS889_II", "BS1308_II", "BS1314_II"), "II", 
                                "III")

cell_freq_stage$stage[cell_freq_stage$grouping == "BS1198_I"] <- "I"
cell_freq_stage$stage[cell_freq_stage$grouping == "BS897_I"] <- "I"
cell_freq_stage$stage[cell_freq_stage$grouping == "BS1318_I"] <- "I"

cell_freq_stage$stage <- factor(cell_freq_stage$stage, levels=c("I","II", "III"))

#add factor levels according to annotation factors
cell_freq$label <- factor(cell_freq$label, levels = c("NK1_FGFBP2", "NK1_CCL4","NK_NFKB1", "NK2",
                                                     "NK3_GZMK", "NK3_CAMK4", "NK3_ENTPD1", "NKprol", 
                                                     "ILC1", "ILC3"))

cell_freq_stage$label <- factor(cell_freq_stage$label, levels = c("NK1_FGFBP2", "NK1_CCL4","NK_NFKB1", "NK2",
                                                                 "NK3_GZMK", "NK3_CAMK4", "NK3_ENTPD1", "NKprol", 
                                                                 "ILC1", "ILC3"))

pdf("to/path/BoxPlot_Histo.pdf", height = 4, width = 8)
cell_freq %>%
  ggplot(aes(label, percent, fill = histo)) +
  geom_boxplot(aes(fill = histo), color = "black", outlier.size = -1, lwd = 0.6) +
  geom_jitter(position=position_jitterdodge(jitter.width = 0.1), aes(color = histo), size = 1) +
  theme_classic() +
  ylim(0,0.5) +
  theme(
    legend.position = "right",
    axis.text.y = element_text(size = 10, colour = "black"),
    axis.text.x = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 13, colour = "black"),
    legend.title = element_text( size = 13),
    legend.text= element_text(size = 13)) +
  labs(x = "",y = "Frequency of NK cells") + 
  scale_fill_manual(values = c("#66CCCC","#FF9999")) + 
  scale_color_manual(values = c("#000000","#000000")) +  
  theme(axis.text.x = element_text(angle = 45,  hjust=1)) +
  geom_pwc(
    aes(group = histo), tip.length = 0.03,
    method = "wilcox_test",
    label.size = 3,
    hide.ns = T,
    bracket.nudge.y =  -0.08, y.position = 0.5, 
    symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns")))
dev.off()

#stage
pdf("to/path/BoxPlot_stage.pdf", height = 4, width = 8)
cell_freq_stage %>%
  ggplot(aes(label, percent, fill = stage)) +
  geom_boxplot(aes(fill = stage), color = "black", outlier.size = -1, lwd = 0.6) +
  geom_jitter(position=position_jitterdodge(jitter.width = 0.1), aes(color = stage), size = 1) +
  theme_classic() +
  ylim(0,0.5) +
  theme(
    legend.position = "right",
    axis.text.y = element_text(size = 10, colour = "black"),
    axis.text.x = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 13, colour = "black"),
    legend.title = element_text(size = 13),
    legend.text= element_text(size = 13)) +
  labs(x = "",y = "Frequency of NK cells") + 
  scale_fill_manual(values = c("#6699FF","#FFCC00", "#CCCCCC")) + 
  scale_color_manual(values = c("#000000","#000000", "#000000")) + 
  theme(axis.text.x = element_text(angle = 45,  hjust=1)) +
  stat_compare_means(
    method = "anova", 
    comparisons = list(c("I", "II"), c("II", "III"), c("I", "III")),
    hide.ns = TRUE,
    tip.length = 0.03,
    label.size = 3) 
dev.off()

# -- Phenotyping Heatmap -- Fig. 2B

#selection of genes assocaited with NK cell phenotype and functions

acR <- c("FCGR3A", "CD69","CD226","ITGB2","NCR1","NCR3","ITGAL",
         "CD244","KLRD1","KLRF1","TNFRSF9","CRTAM", "CD44","KLRK1","CD2","KLRC2","KLRC4",
         "TNFSF10", "NCR2", "TMIGD2","KIR2DL4")

cyto <- c("PRF1","GZMB", "GZMH","TNF","GZMK", "IFNG","GZMA")

TR <- c("CXCR6","ITGAE", "ITGA1", "RGS1","ITGA2")

inR_ICs <- c("KLRG1","CD38","HAVCR2","KIR3DL1","KIR3DL2","KIR2DL3","KIR2DL1",
             "TIGIT", "TOX","ENTPD1","CD96","CTLA4","LAG3","PDCD1","KIR3DL3","KLRC1",  
             "LAYN","KIR3DL1","NT5E")

Proliferation <- c("MKI67", "TOP2A", "NME1")


chem <- c("CCL3", "CCL4",
          "XCL1", "XCL2", "FLT3LG","CCL5","IL32","CXCL13")

cytoR <- c("TGFBR1","TGFBR2","S1PR5","IL12RB1","CX3CR1",
           "CXCR2","CXCR1","IL21R","CCR7","IL18R1","IL12RB2","IL2RB","CXCR3","CXCR4","CCR5",
           "CXCR6")

#select only genes that are expressed more than 10% on the cells within at least one cluster

List_Markers = list.append(acR, cyto, TR, inR_ICs, Proliferation, chem, cytoR)
List_of_Lists= list(acR, cyto, TR, inR_ICs, Proliferation, chem, cytoR)

List_total= intersect(List_Markers, rownames(NK))

GROUP_TO_APPLY_FILTER = "Fine_Annotation"

MIN_PCT_GENE = 10 

perc_exp <- DotPlot(NK, features = List_total, group.by = GROUP_TO_APPLY_FILTER)$data[, c("features.plot", "id", "pct.exp")] 

perc_exp %>%
  filter(pct.exp > MIN_PCT_GENE ) -> Acceptable_Genes 

Acceptable_Genes = Acceptable_Genes$features.plot

#Filter your gene lists
List_total= intersect(List_Markers, Acceptable_Genes)

List_of_Lists=  lapply(List_of_Lists, function(x) intersect(x, Acceptable_Genes))

#now aggregate, exclude ILCs as not needed for question

agg <- AggregateExpression(
  subset(NK, Fine_Annotation =="ILC1"|Fine_Annotation =="ILC3", invert = T),
  assays = "RNA",
  return.seurat = T,
  group.by = "Fine_Annotation")

mat <- agg@assays$RNA$data
colnames <- colnames(mat)
colnames <- gsub('-', '_', colnames)
colnames(mat) <- colnames

mat.scaled= scale(t(mat))

mypalette <- brewer.pal(n = 9, "RdBu")
myBreaks <- seq(2, -2, length.out =9)

#separate matrices
mat1 <- mat.scaled[,intersect(List_of_Lists[[1]], rownames(mat))]
mat1_df <- as.data.frame(mat1)

mat2<- mat.scaled[,intersect(List_of_Lists[[2]], rownames(mat))]
mat2_df <- as.data.frame(mat2)

mat3<- mat.scaled[,intersect(List_of_Lists[[3]], rownames(mat))]
mat3_df <- as.data.frame(mat3)

mat4<- mat.scaled[,intersect(List_of_Lists[[4]], rownames(mat))]
mat4_df <- as.data.frame(mat4)

mat6 <- mat.scaled[,intersect(List_of_Lists[[5]], rownames(mat))]
mat6_df <- as.data.frame(mat6)

mat7<- mat.scaled[,intersect(List_of_Lists[[6]], rownames(mat))]
mat7_df <- as.data.frame(mat7)

mat8<- mat.scaled[,intersect(List_of_Lists[[7]], rownames(mat))]
mat8_df <- as.data.frame(mat8)

#separate heatmaps

h1 <- Heatmap(mat1_df, name = "scaled\nExpression",
              cluster_rows = F,
              show_column_dend = FALSE,
              cluster_row_slices = FALSE,
              row_title_gp = gpar(fontsize = 14),
              cluster_columns = F,
              show_row_dend = FALSE,
              col = colorRamp2(myBreaks, mypalette),
              column_names_gp = gpar(fontsize = 12),
              column_names_side = "bottom",
              column_names_rot = 45,
              show_column_names = T,
              show_row_names = T,
              row_names_side = "left",
              use_raster = F,
              border = T,
              border_gp = gpar(col = "black", lwd = 0.5))

h2 <- Heatmap(mat2_df, name = "scaled\nExpression",
              cluster_rows = F,
              show_column_dend = FALSE,
              cluster_row_slices = FALSE,
              row_title_gp = gpar(fontsize = 14),
              cluster_columns = F,
              show_row_dend = FALSE,
              col = colorRamp2(myBreaks, mypalette),
              column_names_gp = gpar(fontsize = 12),
              column_names_side = "bottom",
              column_names_rot = 45,
              show_column_names = T,
              show_row_names = F,
              row_names_side = "left",
              use_raster = F,
              border = T,
              border_gp = gpar(col = "black", lwd = 0.5))

h3 <- Heatmap(mat3_df, name = "scaled\nExpression",
              cluster_rows = F,
              show_column_dend = FALSE,
              cluster_row_slices = FALSE,
              row_title_gp = gpar(fontsize = 14),
              cluster_columns = T,
              show_row_dend = FALSE,
              col = colorRamp2(myBreaks, mypalette),
              column_names_gp = gpar(fontsize = 12),
              column_names_side = "bottom",
              column_names_rot = 45,
              show_column_names = T,
              show_row_names = F,
              row_names_side = "left",
              use_raster = F,
              border = T,
              border_gp = gpar(col = "black", lwd = 0.5))

h4 <- Heatmap(mat4_df, name = "scaled\nExpression",
              cluster_rows = F,
              show_column_dend = FALSE,
              cluster_row_slices = FALSE,
              row_title_gp = gpar(fontsize = 14),
              cluster_columns = F,
              show_row_dend = FALSE,
              col = colorRamp2(myBreaks, mypalette),
              column_names_gp = gpar(fontsize = 12),
              column_names_side = "bottom",
              column_names_rot = 45,
              show_column_names = T,
              show_row_names = F,
              row_names_side = "left",
              use_raster = F,
              border = T,
              border_gp = gpar(col = "black", lwd = 0.5))

h6 <- Heatmap(mat6_df, name = "scaled\nExpression",
              cluster_rows = F,
              show_column_dend = FALSE,
              cluster_row_slices = FALSE,
              row_title_gp = gpar(fontsize = 14),
              cluster_columns = F,
              show_row_dend = FALSE,
              col = colorRamp2(myBreaks, mypalette),
              column_names_gp = gpar(fontsize = 12),
              column_names_side = "bottom",
              column_names_rot = 45,
              show_column_names = T,
              show_row_names = F,
              row_names_side = "left",
              use_raster = F,
              border = T,
              border_gp = gpar(col = "black", lwd = 0.5))

h7 <- Heatmap(mat7_df, name = "scaled\nExpression",
              cluster_rows = F,
              show_column_dend = FALSE,
              cluster_row_slices = FALSE,
              row_title_gp = gpar(fontsize = 14),
              cluster_columns = F,
              show_row_dend = FALSE,
              col = colorRamp2(myBreaks, mypalette),
              column_names_gp = gpar(fontsize = 12),
              column_names_side = "bottom",
              column_names_rot = 45,
              show_column_names = T,
              show_row_names = F,
              row_names_side = "left",
              use_raster = F,
              border = T,
              border_gp = gpar(col = "black", lwd = 0.5))


h8 <- Heatmap(mat8_df, name = "scaled\nExpression",
              cluster_rows = F,
              show_column_dend = FALSE,
              cluster_row_slices = FALSE,
              row_title_gp = gpar(fontsize = 14),
              cluster_columns = F,
              show_row_dend = FALSE,
              col = colorRamp2(myBreaks, mypalette),
              column_names_gp = gpar(fontsize = 12),
              column_names_side = "bottom",
              column_names_rot = 45,
              show_column_names = T,
              show_row_names = F,
              row_names_side = "left",
              use_raster = F,
              border = T,
              border_gp = gpar(col = "black", lwd = 0.5))


hl = h1+h2+h3+h4+h6+h7+h8

pdf(file = "path/to/Heatmap_phenotyping_10perc.pdf", width =18, height = 2.5)
hl
dev.off()

saveRDS(NK, "path/to/NK_only.RDS")








