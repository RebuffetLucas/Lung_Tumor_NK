library(BayesPrism)
library(Seurat)
library(escape)
library(AUCell)
library(UCell)
library(dplyr)
library(dittoSeq)
library(ggpubr)
library(ggsignif)

# -- 1 - : NDB vs. DCB - Jung et al. 2019, GSE135222:  doi: 10.1038/s41467-019-12159-9 --
# -- 2 - : SD vs. PR vs. CR - Chen et al. 2024, GSE236581:  https://doi.org/10.1016/j.ccell.2024.06.009 --
# -- 3 - : NE vs. E + Pre vs. On - Bassez et al. 2021, (http://biokey.lambrechtslab.org/):  https://doi.org/10.1038/s41591-021-01323-8 --

#  1 : NDB vs. DCB - Jung et al. 2019, GSE135222:  doi: 10.1038/s41467-019-12159-9

#extract NK cells using BayesPrism 
#load image of scDataset as reference prepared in 10_BayesBrism_OverallSurvival.R
load("path/to/Kim_LUAD_sc_tumor.RData")

#load raw bulkMatrix
bulkNSCLC <- readRDS("path/to/Jung_rawMat.RDS")

#row names need to be sample IDs, colnames need to be genes - check + transpose if needed
mt <- as.matrix(mt)
mt <- t(mt)

#visualize outlier genes in bulkRNA
Jung.stat <- plot.bulk.outlier(
  bulk.input=mt, 
  sc.input=sc_lung.filtered.pc, 
  cell.type.labels=cell.type.labels_LUAD,
  species="hs",
  return.raw=TRUE)

#check the concordance of gene expression for different types of genes
plot.bulk.vs.sc(sc.input = sc_lung.filtered.pc,
                bulk.input = mt)


#Contruct a prism object
myPrism_Jung <- new.prism(
  reference=sc_lung.filtered.pc,
  mixture=mt,
  input.type="count.matrix",
  cell.type.labels = cell.type.labels_LUAD,
  cell.state.labels = cell.state.labels_LUAD,
  key="Epithelial cells",
  outlier.cut=0.01,
  outlier.fraction=0.1
)

#RUN PRISM
bp.res_Jung <- run.prism(prism = myPrism_Jung, n.cores=50)


# extract posterior mean of cell type fraction theta
theta_Jung <- get.fraction (bp=bp.res_Jung,
                            which.theta="final",
                            state.or.type="type")

#extract posterior mean of cell type-specific gene expression count matrix Z (NK cells)
NK_Jung <- get.exp(bp=bp.res_Jung,
                   state.or.type="type",
                   cell.name="NK cells")

#transpose to do a seurat object for easy computing
mat_NK_Jung <- t(NK_Jung)

#load metadata of bulkRNA
meta <- read.csv("path/to/meta_GSE135222.csv")

seu <- CreateSeuratObject(
  mat_NK_Jung,
  assay = "RNA",
  names.field = 1,
  names.delim = "_",
  meta.data = meta)

#normalize + scale
seu <- NormalizeData(seu)
seu <- ScaleData(seu)

# - Enrichment analysis - 

#For the bulk --> calculate the gene ranks in order to no which maxrank to use
plotGeneCount(as.matrix(seu@assays$RNA$counts))
cells_rankings <- AUCell_buildRankings(as.matrix(seu@assays$RNA$counts),
                                       plotStats=TRUE)

cells_AUC <- AUCell_calcAUC(list(tumor_NK=top30), cells_rankings)
AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE) 
ceiling(0.15 * nrow(cells_rankings))

#
gene_counts <- rowSums(as.matrix(seu@assays$RNA$counts))

# Plot the histogram

p1 <- ggplot(data.frame(Counts = gene_counts), aes(x = Counts)) +
  geom_histogram(fill = "skyblue", color = "black", alpha = 0.7) +
  scale_x_log10() +  # Log scale to handle wide range of counts
  ggtitle("Distribution of Counts Per Feature") +
  xlab("Counts in log scale") +
  ylab("Number of Features (Genes)") +
  theme_minimal()

p2 <- ggplot(data.frame(Counts = gene_counts), aes(x = Counts)) +
  geom_histogram(binwidth = 500, fill = "skyblue", color = "black", alpha = 0.7) +
  #scale_x_log10() +  # Log scale to handle wide range of counts
  ggtitle("Distribution of Counts Per Feature") +
  xlab("Counts") +
  ylab("Number of Features (Genes)") +
  theme_minimal()

p1+p2
#select maxRank of 2000

#load signature to test
tNK <- read.csv("/path/to/taNK_markers.csv")
tNK <- tNK[order(tNK$avg_log2FC, decreasing = T),]
top30 <- tNK$X[1:30] #top30 genes

#enrichment
ES <- enrichIt(obj = as.matrix(seu@assays$RNA$counts), 
               gene.sets = list(tumor_NK_30=top30), 
               groups = 1000, cores = 12, maxRank = 2000,
               min.size = 5, method = "UCell")

seu <- AddMetaData(seu,ES)
ES <- data.frame(seu[[]])

ES$clinical.benefit <- factor(ES$clinical.benefit, levels = c("NDB", "DCB"))

#again dittoseq can be used for easy plotting
counts <- ES %>% dplyr::select(tumor_NK_30) %>% as.matrix()
counts <- t(counts)
metadata <- ES %>% dplyr::select(-tumor_NK_30)
bulkSCE <- importDittoBulk(
  x = list(counts = counts),
  metadata = data.frame(metadata))

Colors <- c("#6699FF", "#FF6666")

pdf("/path/to/UCell_Box.pdf", width = 2, height = 2.5)
dittoBoxPlot(bulkSCE, "tumor_NK_30", group.by = "clinical.benefit", color.panel = Colors, 
             boxplot.width = 0.5, 
             boxplot.color = "black", boxplot.lineweight = 0.3, min = 0, max = 0.08,
             ylab = "tumor-associated NK \nsignature score", jitter.size = 1.2) +
  geom_pwc(group.by = "x.var", hide.ns = TRUE, method = "wilcox_test", step.increase = 0.2, label = "p.signif",
           symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), 
                              symbols = c("****", "***", "**", "*", "ns"))) &
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.text.y = element_text(size = 9, color = "black"),
        axis.text.x = element_text(size = 13, color = "black", angle = 45),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 11, color = "black", face = "bold"),
        legend.position = "none") 
dev.off()

# -- 2 - : SD vs. PR vs. CR - Chen et al. 2024, GSE236581:  https://doi.org/10.1016/j.ccell.2024.06.009 --

#load the dataset, extract NK cells

CRC <- Read10X(data.dir = "path/to/Chen_GSE236581/data/")
meta <- read.csv("path/to/Chen_GSE236581/CRC_meta.csv")
rownames(meta) <- meta$X

# Initialize the Seurat object with the raw (non-normalized data).
CRC <- CreateSeuratObject(counts = CRC, meta.data =  meta, min.cells = 3, min.features = 200)

#subset tumor tissue and innate lymphoid cells to get NK cells
CRC_subset <- subset(CRC, Tissue == "Tumor" & 
                       MajorCellType == "ILC")

#check QC
VlnPlot(CRC_subset, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3)
CRC_subset <- subset(CRC_subset,nCount_RNA<10000)

#add response information according to paper supplementary information
CRC_subset$Response <- ifelse(CRC_subset$Patient %in% c("P24","P17", "P14", 
                                                        "P15", "P05", "P03"), "PR", "CR")

CRC_subset$Response[CRC_subset$Patient %in% c("P02", "P16", "P25")] <- "SD"

#standart Seurat workflow
CRC_subset <- NormalizeData(CRC_subset)
CRC_subset <- FindVariableFeatures(CRC_subset, nfeatures = 1000)
CRC_subset <- ScaleData(CRC_subset)
CRC_subset <- RunPCA(CRC_subset)
#ElbowPlot(CRC_subset)

CRC_subset <- FindNeighbors(CRC_subset, dims = 1:15)
CRC_subset <- FindClusters(CRC_subset, resolution = seq(0, 1, 0.5))
CRC_subset <- RunUMAP(CRC_subset, dims = 1:15)

#DimPlot(CRC_subset, group.by = "RNA_snn_res.0.5")
CRC_subset <- SetIdent(CRC_subset, value = CRC_subset$RNA_snn_res.0.5)

CRC_subset$RNA_snn_res.0.5 <- as.numeric(CRC_subset$RNA_snn_res.0.5)
Markers_CRC <- FindAllMarkers(CRC_subset, assay = "RNA", min.pct = 0.25)

Markers_CRC  %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) -> top20

#Annotate

CRC_subset$Annotation <- "NK1_FCGR3A"
CRC_subset$Annotation[CRC_subset$RNA_snn_res.0.5 == 1] <- "NK2_GZMK"
CRC_subset$Annotation[CRC_subset$RNA_snn_res.0.5 == 3] <- "NK2"
CRC_subset$Annotation[CRC_subset$RNA_snn_res.0.5 == 4] <- "NK3_ITGAE"
CRC_subset$Annotation[CRC_subset$RNA_snn_res.0.5 == 5] <- "ILC3"
CRC_subset$Annotation[CRC_subset$RNA_snn_res.0.5 == 6] <- "NKstressed"
CRC_subset$Annotation[CRC_subset$RNA_snn_res.0.5 == 7] <- "NK1_CX3CR1"
CRC_subset$Annotation[CRC_subset$RNA_snn_res.0.5 == 8] <- "NK3_CXCR6"
CRC_subset$Annotation[CRC_subset$RNA_snn_res.0.5 == 9] <- "NKprol"

CRC_subset <- SetIdent(CRC_subset, value = CRC_subset$Annotation)
CRC_subset$Annotation <- factor(CRC_subset$Annotation, levels = c("NK1_FCGR3A","NK1_CX3CR1",
                                                                "NK2","NK2_GZMK","NK3_CXCR6",
                                                                "NK3_ITGAE","NKprol",
                                                                "NKstressed","ILC3"))

Colors <- c("#E15759","#FF9D9A",  "#59A14F", "#8CD17D", "#9966CC", "#CC99FF","#499894", "#86BCB6", "#6495ED")


pdf("path/to/Chen_GSE236581/Dimplot.pdf", width = 4.5, height = 3) #- Fig. S6C
DimPlot(CRC_subset, group.by = "Annotation", cols = Colors) + NoAxes()
dev.off()

pdf(file = "path/to/Chen_GSE236581/DotPlot_Annotated.pdf", width = 7, height = 2.5) #- Fig. S6D
DotPlot(CRC_subset, group.by = "Annotation",
        features =  c("NKG7","FCGR3A", "CX3CR1","NCAM1","SELL","XCL1", "CD44",
                      "GZMK", "CXCR3", "CD52", "IL32",
                      "ITGA1", "ITGAE", "ENTPD1", "GZMA","CXCR6","MKI67", "TOP2A",
                      "DNAJB1", "HSPA1B", "RORC", "IL1R1"),
        assay= "RNA", cols = c("RdBu")) + RotatedAxis() +
  theme(axis.text.y = element_text(size = 9),
        axis.text.x = element_text(size = 9, face = "bold"),
        legend.title = element_text(size=9),
        legend.text = element_text(size=8),
        legend.key.size = unit(0.3, "cm"),
        axis.title = element_blank())
dev.off()

#enrichment
agg_NK <- AggregateExpression(CRC_subset, assays = "RNA",return.seurat = T, group.by = c("Response", "Patient"))
agg_NK_noseu <- AggregateExpression(CRC_subset, assays = "RNA",return.seurat = F, group.by = c("Response", "Patient"))

ES_UCELL <- enrichIt(obj = as.matrix(agg_NK_noseu[["RNA"]]), 
                     gene.sets =  list(tumor_NK_30=top30), 
                     groups = 1000, cores = 12, 
                     min.size = 5, method = "UCell")



agg_NK <- AddMetaData(agg_NK, ES_UCELL)
ES_UCell_df <- data.frame(agg_NK[[]])

ES_UCell_df$Response <- factor(ES_UCell_df$Response, levels = c("SD", "PR", "CR"))

counts <- ES_UCell_df %>% dplyr::select(tumor_NK_30) %>% as.matrix()
counts <- t(counts)
metadata <- ES_UCell_df %>% dplyr::select(-tumor_NK_30)

bulkSCE<- importDittoBulk(
  x = list(counts = counts),
  metadata = data.frame(metadata))

Colors <- c("#6699FF", "#FFCC33", "#FF6666")

pdf("path/to/UCell_Box_top30.pdf", width = 2, height = 2.5)  #- Fig. S6E
dittoBoxPlot(bulkSCE, "tumor_NK_30", group.by = "Response", color.panel = Colors, 
                    boxplot.width = 0.5, 
                    boxplot.color = "black",boxplot.lineweight = 0.3, min = 0, max = 0.2,
                    ylab = "tumor-associated NK \nsignature score", jitter.size = 1.2) +
  geom_pwc(group.by="x.var",label="p.adj.signif",hide.ns=T, method = "dunn_test", step.increase = 0.2) &
  theme(plot.title = element_text(size = 12,face = "bold", hjust = 0.5),
        axis.text.y = element_text(size = 9,color = "black"),
        axis.text.x = element_text(size = 13,  color = "black", angle = 45),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 11, color = "black", face = "bold"),
        legend.position = "none")
dev.off()


# -- 3 - : NE vs. E + Pre vs. On - Bassez et al. 2021, (http://biokey.lambrechtslab.org/):  https://doi.org/10.1038/s41591-021-01323-8 --

C1 <- readRDS("path/to/1863-counts_cells_cohort1.rds")
C2 <- readRDS("path/to/1867-counts_cells_cohort2.rds")
Met2 <- read.csv("path/to/C2.csv")
Met1 <-  read.csv("path/to/C1.csv")

S1 <- CreateSeuratObject(C1,assay = "RNA",  meta.data = Met1)
S2 <- CreateSeuratObject(C2,assay = "RNA",  meta.data = Met2)

merge <- merge(S1, S2, merge.data = T)

#seurat workflow
merge <- NormalizeData(merge)
merge <- FindVariableFeatures(merge, selection.method = "vst", nfeatures = 1500)
merge <- ScaleData(merge)
merge <- RunPCA(merge, features = VariableFeatures(object = merge))

#ElbowPlot(merge)
merge <- FindNeighbors(merge, dims = 1:15)
merge <- FindClusters(merge, resolution = seq(0.1, 1, 0.1))
merge <- RunUMAP(merge, dims = 1:15)

#extract NK cells
#check for expression of main marker
FeaturePlot(merge, c("CD3D", "CD3G", "CD4", "CD8A", "NKG7"))

#only take lymphocyte marker and remove other celltypes from annotation

sub <- subset(merge, RNA_snn_res.0.1 == 13|
                RNA_snn_res.0.1 == 0|
                RNA_snn_res.0.8 == 28)

sub <- subset(sub, cellType == "B_cell" |
                  cellType == "Cancer_cell" |
                  cellType == "Endothelial_cell" |
                  cellType == "Fibroblast" |
                  cellType == "Myeloid_cell" |
                  cellType == "pDC"|
                  cellType == "Mast_cell", 
                  invert = T)

#remove T cells
RidgePlot(sub, c("CD3D", "CD3G", "CD4", "CD8A"))

expr <- FetchData(sub, vars = c("CD3D", "CD3G", "CD8A"))
cells_keep <- rownames(expr)[expr$CD3D <= 1 & expr$CD3G <= 0.5 & expr$CD8A <= 0.5]
sub_NK <- subset(sub, cells = cells_keep)

#reanalyse
sub_NK[["percent.mt"]] <- PercentageFeatureSet(sub_NK, pattern = "^MT-")
sub_NK[["percent.rb"]] <- PercentageFeatureSet(sub_NK, pattern = "^RP[S|L]")

sub_NK <- NormalizeData(sub_NK)
sub_NK <- FindVariableFeatures(sub_NK, selection.method = "vst", nfeatures = 1000)

sub_NK <- ScaleData(sub_NK, vars.to.regress = c("percent.mt", "percent.rb"))
sub_NK <- RunPCA(sub_NK, features = VariableFeatures(object = sub_NK))

#ElbowPlot(sub_NK)
sub_NK <- FindNeighbors(sub_NK, dims = 1:10)
sub_NK <- FindClusters(sub_NK, resolution = seq(0.1, 0.5, 0.1))
sub_NK <- RunUMAP(sub_NK, dims = 1:10)

#still some non-NK cells, further clean until pure
NK_bassez <- subset(sub_NK, 
                    RNA_snn_res.0.5 == 0 |
                    RNA_snn_res.0.5 == 1 |
                    RNA_snn_res.0.5 == 3 |
                    RNA_snn_res.0.5 == 5, 
                    invert = T)

NK_bassez <- NormalizeData(NK_bassez)
NK_bassez <- FindVariableFeatures(NK_bassez, selection.method = "vst", nfeatures = 1000)

NK_bassez <- ScaleData(NK_bassez, vars.to.regress = c("percent.mt", "percent.rb"))
NK_bassez <- RunPCA(NK_bassez, features = VariableFeatures(object = NK_bassez))

#ElbowPlot(NK_bassez)
NK_bassez <- FindNeighbors(NK_bassez, dims = 1:10)
NK_bassez <- FindClusters(NK_bassez, resolution = seq(0.1, 0.5, 0.1))
NK_bassez <- RunUMAP(NK_bassez, dims = 1:10)

Markers_NK <- FindAllMarkers(NK_bassez, assay = "RNA", min.pct = 0.25, test.use = "MAST")

Markers_NK  %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) -> top20

#Annotate
DimPlot(NK_bassez, group.by = "RNA_snn_res.0.5")
NK_bassez$RNA_snn_res.0.5 <- as.numeric(NK_bassez$RNA_snn_res.0.5)

NK_bassez$Annotation <- "NK2"

NK_bassez$Annotation[NK_bassez$RNA_snn_res.0.5 == 2] <- "NK3_ITGAE"
NK_bassez$Annotation[NK_bassez$RNA_snn_res.0.5 == 3] <- "NK1"
NK_bassez$Annotation[NK_bassez$RNA_snn_res.0.5 == 4] <- "ILC1"
NK_bassez$Annotation[NK_bassez$RNA_snn_res.0.5 == 5] <- "NKprol"
NK_bassez$Annotation[NK_bassez$RNA_snn_res.0.5 == 6] <- "NK_ISGs"
NK_bassez$Annotation[NK_bassez$RNA_snn_res.0.5 == 7] <- "ILC3"
NK_bassez$Annotation[NK_bassez$RNA_snn_res.0.5 == 8] <- "NK1"


NK_bassez$Annotation <- factor(NK_bassez$Annotation, levels = c("NK1",
                                                                "NK2",
                                                                "NK3_ITGAE",
                                                                "NK_ISGs",
                                                                "NKprol",
                                                                "ILC1",
                                                                "ILC3"))

NK_bassez <- SetIdent(NK_bassez, value = NK_bassez$Annotation)

Colors <- c("#E15759", "#59A14F", "#8CD17D", "#9966CC", "#CC99FF","#499894", "#86BCB6", "#6495ED")

pdf(file = "path/to/DimPlot.pdf",width = 4, height = 2.5)  #Fig. S6F
DimPlot(NK_bassez, cols = Colors) & NoAxes()
dev.off()

pdf(file = "path/to/Dotplot.pdf",width = 6, height = 2.5) #Fig. S6F
DotPlot(NK_bassez, features = c("FCGR3A", "CX3CR1", "PRF1", "NKG7", "XCL1", "XCL2", "SELL", "CD52", "IL32",
                                "ITGAE", "CCL5","ITGA1","ENTPD1", "IFIT1", "ISG15", "IRF7", "MKI67", "TOP2A", 
                                "RORA", "CD6", "KIT", "RORC"),
        assay= "RNA", cols = c("RdBu")) + RotatedAxis()  + theme(axis.text.y = element_text(size = 9),
                                                                 axis.text.x = element_text(size = 9, face = "bold"),
                                                                 axis.title = element_blank(),
                                                                 legend.title = element_text(size=9),
                                                                 legend.text = element_text(size=8),
                                                                 legend.key.size = unit(0.3, "cm"))
dev.off()

#Enrichment - time= timepoinst, Exp = Expander
agg_NK_time <- AggregateExpression(NK_bassez, assays = "RNA",return.seurat = T, group.by = c("timepoint", "patient_id"))
agg_NK_Exp <- AggregateExpression(NK_bassez, assays = "RNA",return.seurat = T, group.by = c("expansion", "patient_id"))

agg_NK_time_noseu <- AggregateExpression(NK_bassez, assays = "RNA",return.seurat = F, group.by = c("timepoint", "patient_id"))
agg_NK_Exp_noseu <- AggregateExpression(NK_bassez, assays = "RNA",return.seurat = F, group.by = c("expansion", "patient_id"))

ES_time <- enrichIt(obj = as.matrix(agg_NK_time_noseu[["RNA"]]), 
                    gene.sets = list(tumor_NK_30=top30), 
                    groups = 1000, cores = 12, 
                    min.size = 5, method = "UCell")


agg_NK_time <- AddMetaData(agg_NK_time, ES_time)
ES_time_df <- data.frame(agg_NK_time[[]])

ES_time_df$timepoint <- factor(ES_time_df$timepoint, levels = c("Pre", "On"))

ES_Exp <- enrichIt(obj = as.matrix(agg_NK_Exp_noseu[["RNA"]]), 
                   gene.sets = list(tumor_NK_30=top30), 
                   groups = 1000, cores = 12, 
                   min.size = 5, method = "UCell")


agg_NK_Exp <- AddMetaData(agg_NK_Exp, ES_Exp)
ES_Exp_df <- data.frame(agg_NK_Exp[[]])

ES_Exp_df$expansion <- factor(ES_Exp_df$expansion, levels = c("NE", "E"))

#TimePoint
counts <- ES_time_df %>% dplyr::select(tumor_NK_30) %>% as.matrix()
counts <- t(counts)
metadata <- ES_time_df %>% dplyr::select(-tumor_NK_30)

bulkSCE_Time <- importDittoBulk(
  x = list(counts = counts),
  metadata = data.frame(metadata))

Colors <- c("#6699FF", "#FFCC33")

pdf("path/to/Box_Timepoint_top30.pdf", width = 2, height = 2.5) #Fig. S6G
dittoBoxPlot(bulkSCE_Time, "tumor_NK_30", group.by = "timepoint", color.panel = Colors, 
                   boxplot.width = 0.5, 
                   boxplot.color = "black",boxplot.lineweight = 0.3, min = 0, max = 0.13,
                   ylab = "tumor-associated NK \nsignature score", jitter.size = 1.2) +
  geom_pwc(group.by="x.var",label="p.adj", hide.ns = T,method = "wilcox_test", step.increase = 0.2) &
  theme(plot.title = element_text(size = 13,  color = "black", face = "bold", hjust = 0.5),
        axis.text.y = element_text(size = 9,color = "black"),
        axis.text.x = element_text(size = 13,  color = "black", angle = 45),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 11, color = "black", face = "bold"),
        legend.position = "none") &ggtitle("top30 genes")
dev.off()


#Expansion
ES_Exp_df <- na.omit(ES_Exp_df)
counts <- ES_Exp_df %>% dplyr::select(tumor_NK_30) %>% as.matrix()
counts <- t(counts)
metadata <- ES_Exp_df %>% dplyr::select(-tumor_NK_30)
bulkSCE_Exp <- importDittoBulk(
  x = list(counts = counts),
  metadata = data.frame(metadata))

Colors <- c("#499894","#FF6666")


pdf("path/to/Box_Expander_top30.pdf", width = 2, height = 2.5) #Fig. S6G
dittoBoxPlot(bulkSCE_Exp, "tumor_NK_30", group.by = "expansion", color.panel = Colors, 
                   boxplot.width = 0.5, 
                   boxplot.color = "black",boxplot.lineweight = 0.3, min = 0, max = 0.13,
                   ylab = "tumor NK \nsignature score", jitter.size = 1.2) +
  geom_pwc(group.by="x.var", hide.ns = T,method = "wilcox_test", step.increase = 0.4,
           label = "p.signif",
           symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))) &
  theme(plot.title = element_text(size = 13,  color = "black", face = "bold", hjust = 0.5),
        axis.text.y = element_text(size = 9,color = "black"),
        axis.text.x = element_text(size = 13,  color = "black", angle = 45),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 11, color = "black", face = "bold"),
        legend.position = "none") &ggtitle("top 30 genes")
dev.off()




























