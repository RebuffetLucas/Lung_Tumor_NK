library(Seurat)
library(SeuratDisk)
library(ggpubr)
library(dittoSeq)
library(ggrepel)
library(dplyr)
library(escape) #v1.12.0
set.seed(123456)

# -- Preparation of non-tumor lung NK cell datasets --

# - 1. Bischoff et. al. - DOI: https://doi.org/10.1038/s41388-021-02054-3 - Download: https://codeocean.com/capsule/8321305/tree/v1 -

#only normal samples were downloaded (p...n)

patients <- list()
p <- list("18","19", "27", "28", "29", "30", "31", "32", "33", "34")
counts <- list()

for (i in p){
  counts[[i]] <- Read10X(paste("~/healthy_lung/cellranger/normal/p0",i,"n/filtered_feature_bc_matrix", sep = ""))
  patients[[i]] <- CreateSeuratObject(
    counts = counts[[i]],
    min.cells = 3, min.features = 200,
    assay = "RNA",
    project = paste("p0", i, "n", sep = ""))
}

#addd metadata - information obtained from online platform
patients[[1]]$sex <- "female"
patients[[2]]$sex <- "male"
patients[[3]]$sex <- "male"
patients[[4]]$sex <- "male"
patients[[5]]$sex <- "female"
patients[[6]]$sex <- "male"
patients[[7]]$sex <- "male"
patients[[8]]$sex <- "female"
patients[[9]]$sex <- "female"
patients[[10]]$sex <- "male"

#QC
integrated <- merge(patients[[1]], c(patients[2:10]), 
                add.cell.ids = names(patients), merge.data = T)

integrated[["percent.rb"]] <- PercentageFeatureSet(integrated, pattern = "^RP[S|L]")
integrated[["percent.mt"]] <- PercentageFeatureSet(integrated, pattern = "^MT-")

QC.features <- c("nFeature_RNA", "nCount_RNA","percent.mt",
                 "percent.rb")

VlnPlot(integrated, features = QC.features, group.by = "orig.ident",
        ncol = 3, pt.size = 0.1)


#SUBSET
integrated <- subset(
  x = integrated,
  subset = nCount_RNA < 50000 &
    nCount_RNA > 1000 &
    nFeature_RNA > 600 &
    percent.mt < 25 &
    percent.rb > 0.05
)

#Normalization + Integration
list <- SplitObject(integrated, split.by = "orig.ident")

list <- lapply(
  X = list,
  FUN = SCTransform,
  assay = "RNA",
  verbose = T,
  vst.flavor = "glmGamPoi",
  vars.to.regress = c("percent.mt", "percent.rb"),
  return.only.var.genes = F
)

features <- SelectIntegrationFeatures(list)
gc()
list <- PrepSCTIntegration(object.list = list, anchor.features = features)
gc()
anchors <- FindIntegrationAnchors(object.list = list, normalization.method = "SCT",
                                  anchor.features = features)
gc()
integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

DefaultAssay(integrated) <- "integrated"
integrated  <- RunPCA(integrated, npcs = 75, verbose = T)
ElbowPlot(integrated, ndims = 75) + geom_hline(yintercept=2, linetype="dashed", color = "red") #PC30

integrated <- FindNeighbors(integrated, dims = 1:30)
integrated <- FindClusters(integrated, graph.name = "integrated_snn", resolution = seq(0.1, 1, by=0.1))
integrated <- RunUMAP(integrated, dims = 1:30, verbose = T, min.dist = 0.1)

integrated <- Seurat::SetIdent(integrated, value = integrated$integrated_snn_res.0.8)

DefaultAssay(integrated) <- "SCT"

#check NK and T cell genes
FeaturePlot(integrated, c("NCAM1","NKG7", "GNLY", "ITGAE", "ITGA1", "ITGA2", "CD3E", "CD3D", "CD3G", 
                          "CD8A"), min.cutoff = 0, ncol = 3) & NoAxes() & NoLegend()

#only take clusters with "NKG7", "GNLY" expression

#only lymphocytes taken
subset <- subset(integrated, 
                 integrated_snn_res.0.5 == 4 |
                 integrated_snn_res.0.5 == 8 |
                 integrated_snn_res.0.5 == 22)

# was performed three times until cleaned based on NK and T cells genes, only NK cells derived, characterized by positive expression of
#NKG7, NCAM1 and GNLY and negative expression for CD3D and CD3G

saveRDS(integrated, "path/to/Bischoff_NK_H.RDS")

# - 2. Sikkema et. al. - DOI: https://doi.org/10.1038/s41591-023-02327-2 
# - Download of Core Dataset: https://cellxgene.cziscience.com/collections/6f6d381a-7701-4781-935c-db10d30de293 -

#convert h5ad
Convert("~/Downloads/db3a3799-3c3f-44db-96a5-81da79a5f3e3.h5ad", dest = "h5seurat", overwrite = TRUE)
seurat_obj <- LoadH5Seurat("path/to/db3a3799-3c3f-44db-96a5-81da79a5f3e3.h5seurat")
saveRDS(seurat_obj, file = "Seurat_lungAtlas.rds")

#select only normal tisseu and NK cells 
lungAtlas <- readRDS("path/to/Seurat_lungAtlas.rds")

NK_lung <- subset(lungAtlas, 
                  disease == "normal" &
                  cell_type == "natural killer cell" & 
                  tissue  == "lung parenchyma")

saveRDS(NK_lung, "path/to/Sikkema_NK_H.RDS")


# -- Integration with tumor NK cells --

Sikkema <- readRDS("path/to/Sikkema_NK_H.RDS")
Bischoff <- readRDS("path/to/Bischoff_NK_H.RDS")
NKT <- readRDS("path/to/NK_only.RDS")

#add important information
Sikkema$Origin <- "nontumor_lung"
Sikkema$Paper <- "Sikkema_HumanLungAtlas_2023"
Sikkema$study <- "Sikkema_2023"

Bischoff$study <- "Bischoff_2021"
Bischoff$Origin <- "nontumor_lung"
Bischoff$Paper <- "Bischoff_2021"

NKT$study <- "Zippelius"
NKT$Origin <- "tumor_lung"
NKT$Paper <- NA

#split by donor or/and study
NK_T_list <- SplitObject(NKT, split.by = "orig.ident")
Bischoff_list <- SplitObject(Bischoff, split.by = "orig.ident")
Sikkema_list <- SplitObject(Sikkema, split.by = "study")

#in H2_list, studies with low amount of cells --> remove
Sikkema_list$Misharin_Budinger_2018 <- NULL
Sikkema_list$Misharin_2021 <- NULL
Sikkema_list$Teichmann_Meyer_2019 <- NULL

Sikkema_list$orig.ident <- Sikkema_list$study

#merge all
combined_list <- c(NK_T_list, Bischoff_list, Sikkema_list)
sample_names <- names(combined_list)

integrated <- merge(combined_list[[1]],combined_list[2:length(combined_list)],
  add.cell.ids = sample_names,
  merge.data = TRUE)

#Normalize and integrate
DefaultAssay(integrated) <- "RNA"
integrated[["RNA"]] <- split(integrated[["RNA"]], f = integrated$orig.ident)
integrated <- NormalizeData(integrated)
integrated <- FindVariableFeatures(integrated, nfeatures = 1000)

features <- VariableFeatures(integrated)

load("/scicore/home/zippeliu/serger0000/Multiome/analysis/data/exclude.gene.misc.human.v4.RData")
f.feat <- !(features %in% all.gene.ignore.df[["seu.id"]]) &
  !(grepl("^(RP[LS]|Rp[ls])", features, perl=T)) &
  !(grepl("^MT-|^MTRN|^MTAT|^MTND|^MRP",features,perl=T)) &
  !(grepl("^RP[0-9]+-",features,perl=T)) &
  !(features %in% c('XIST')) 

features.flt <-features[f.feat]

VariableFeatures(integrated) <- features.flt

integrated <- ScaleData(integrated, features = rownames(integrated), vars.to.regress = c("percent.mt"))
integrated <- RunPCA(integrated)

integrated <- IntegrateLayers(object = integrated, method = CCAIntegration, k.weight = 55,
                          orig.reduction = "pca", new.reduction = "integrated.cca",
                          verbose = FALSE)


integrated[["RNA"]] <- JoinLayers(integrated[["RNA"]])

#clustering

#ElbowPlot(integrated)
integrated <- FindNeighbors(integrated, dims = 1:15, reduction = "integrated.cca")
integrated <- RunUMAP(integrated, dims = 1:15, reduction.name = "umap",reduction = "integrated.cca")
integrated <- FindClusters(integrated,resolution = seq(0, 1.5, 0.1))

plt <- list()
y <- seq(0.1,1.5,by=0.1)

for (i in 1:15){
  plt[[i]] <- DimPlot(integrated, reduction = "umap", label = T, label.size = 6,
                      group.by = paste("RNA_snn_res.", y[[i]], sep = "")) + ggtitle(paste("resolution", y[[i]], sep = " ")) &
    NoLegend() & NoAxes()
}

pdf(file = "path/to/res.pdf",
    height = 9, width = 15)
ggarrange(plotlist = plt, nrow = 3,ncol = 5)
dev.off()

#res 0.2
integrated <- SetIdent(integrated, value = integrated$RNA_snn_res.0.2)

#Markers
Markers <- FindAllMarkers(integrated, min.pct = 0.25, assay = "RNA", test.use = "MAST")

Markers  %>%
  group_by(cluster) %>%
  top_n(n = 30, wt = avg_log2FC) -> top30

write.csv(Markers, "path/to/Markers.csv")
write.csv(top30, "path/to/Markers_top30.csv")

#Annotation based on Marker genes - cluster annotations depend on marker expression
integrated$Annotation <- "NK1"
integrated$Annotation[integrated$RNA_snn_res.0.2 == 2] <- "taNK3"
integrated$Annotation[integrated$RNA_snn_res.0.2 == 3] <- "NK3"
integrated$Annotation[integrated$RNA_snn_res.0.2 == 4] <- "NK2"
integrated$Annotation[integrated$RNA_snn_res.0.2 == 5] <- "ILC3"
integrated$Annotation[integrated$RNA_snn_res.0.2 == 6] <- "ILC1"
integrated$Annotation[integrated$RNA_snn_res.0.2 == 7] <- "NKprol"

#make factors from 
integrated$Annotation <- factor(integrated$Annotation, levels = c("NK1", "NK2", "NK3",
                                                          "taNK3", "NKprol", "ILC1","ILC3"))
integrated <- SetIdent(integrated, value = integrated$Annotation)

#Plot Origin and Study
Colors <- c("#E15759","#86BCB6","#59A14F","#6699FF", "#CC99FF","#999999", "#CCCCCC")

#Origin - Fig. 2C
pdf("path/to/DimPlot_split_Origin.pdf", width = 7.5, height = 4)
DimPlot(integrated, reduction = "umap", label = F, group.by = "Annotation",cols = Colors, split.by = "Origin",
        label.size = 4,label.box = T, repel = T, raster = F) & NoAxes()
dev.off()

#Paper - Fig. S3A
pdf("path/to/DimPlot_split_Paper.pdf", width = 5.5, height = 6)
DimPlot(integrated, reduction = "umap", label = F, group.by = "Annotation",cols = Colors, split.by = "Paper",
        label.size = 4,label.box = T, repel = T, raster = F, ncol = 2) & NoAxes() & NoLegend()
dev.off()


#DotPlot - Fig. S3A
pdf(file = "path/to/DotPlot_Anno.pdf", width =6.5, height = 2.5)
DotPlot(integrated, group.by = "Annotation",
        features =  c("NKG7","FCGR3A", "FGFBP2","NCAM1","IL7R", "XCL1", "GZMH","IL32","CD3E", "CD52","GZMK",
                      "ITGAE", "ITGA1", "ENTPD1", "CAMK4","MKI67", "TOP2A","RORA", "CD5", "CD6", "IL1R1", "AHR", "RORC"),
        assay= "RNA", cols = c("RdBu")) + RotatedAxis() +
  theme(axis.text.y = element_text(size = 9),
        axis.text.x = element_text(size = 9, face = "bold"),
        legend.title = element_text(size=9),
        legend.text = element_text(size=8),
        legend.key.size = unit(0.3, "cm"),
        axis.title = element_blank())
dev.off()

# --- Signature Scoring - Fig. S3C ---
signatures <- readRDS("path/to/signatures_List.RDS") #List of signatures to test obtained from indicated papers
#13-gene NK signature, NK1, NK2, NK3, ILCvsNK, ILC1, ILC2, ILC3

agg <- AggregateExpression(integrated, assays = "RNA",return.seurat = T,  slot = "scale.data",
                           group.by = c("Paper", "Annotation"))

agg_noseu <- AggregateExpression(integrated, assays = "RNA",return.seurat = F, slot = "scale.data",
                                 group.by = c("Paper", "Annotation"))


agg@meta.data
agg$cluster <- gsub('-', '_', agg$Annotation)

agg$cluster <- factor(agg$cluster, levels = levels(integrated$Annotation))

#enrich signatures with UCELL
ES <- enrichIt(obj = as.matrix(agg_noseu[["RNA"]]), 
               gene.sets = signatures, 
               groups = 1000, cores = 12, 
               min.size = 5, method = "UCell")

agg <- AddMetaData(agg, ES)
ES <- data.frame(agg[[]])

#fake a dittobulk object for easy plotting
Main_names <- names(signatures)[1:8] # select only the ones you want to plot (NK1-3, ILC sigs)

counts <- ES %>% select(-Annotation, -Paper, -orig.ident, -cluster) %>% as.matrix()
counts <- t(counts)
metadata <- ES %>% select(Paper, cluster)

bulkSCE <- importDittoBulk(
  x = list(counts = counts),
  metadata = data.frame(metadata))

max <- c(1.5, 1.8, 1.8, 1.3,1, 1, 0.6, 1.5) #adjusted to y axis

y <- c("13-gene NK signature score","NK1 signature score", "NK2 signature score", "NK3 signature score",
       "ILC signature score", "ILC1 signature score","ILC2 signature" ,"ILC3 signature")
plt <- list()

for (i in 1:length(Main_names)){
  plt[[i]] <- dittoBoxPlot(bulkSCE, Main_names[[i]], group.by = "cluster", color.panel = Colors, 
                           boxplot.width = 0.5,
                           boxplot.color = "black",boxplot.lineweight = 0.3, min = 0, #max = max[[i]],
                           ylab = y[[i]], jitter.size = 1.2, main = y[[i]]) &
    geom_pwc(group.by="x.var",label="p.adj.signif",hide.ns=TRUE, method = "dunn_test", step.increase = 0.2) &
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 9),
          axis.text.y = element_text(size = 7,color = "black"),
          axis.text.x = element_text(size = 7,  color = "black", angle = 45, face = "bold"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 7, color = "black"),
          legend.position = "none")
}

pdf("path/to/Boxplot_sigscoring.pdf", 
    width =8, height = 4)
ggarrange(plotlist = plt, nrow = 2,ncol = 4)
dev.off()

#Abundacy BoxPlot - Fig. S3D

#add all donor_ids into orig.ident, for the Bischoff and Zippelius study, already done
integrated$orig.ident[is.na(integrated$orig.ident) | integrated$orig.ident == ""] <- 
  integrated$donor_id[is.na(integrated$orig.ident) | integrated$orig.ident == ""]

#retrieve table with Annotation per donor/orig.ident
table <- dittoBarPlot(integrated, var = "Annotation", group.by = "orig.ident", retain.factor.levels = F,data.out=T)
cell_freq <- table$data
cell_freq$Origin <- ifelse(cell_freq$grouping %in% c("BS1175", "BS1198", "BS1308", "BS1314",
                                                     "BS1318", "BS1319", "BS1322", "BS824", "BS840", "BS889", "BS897"), "tumor lung", "non-tumor lung")
cell_freq$Origin <- factor(cell_freq$Origin, levels=c("tumor lung", "non-tumor lung"))
cell_freq$label <- as.character(cell_freq$label)
cell_freq$cluster <- cell_freq$label
cell_freq$label <- factor(cell_freq$label, levels=levels(integrated$Annotation))
cell_freq$cluster <- factor(cell_freq$cluster, levels = c("NK1","NK2","NK3","taNK3","NKprol","ILC1","ILC3" ))

pdf("path/to/Frequency_NKs.pdf", height = 3, width = 4)
cell_freq %>%
  ggplot(aes(cluster, percent, color = Origin)) +
  geom_boxplot(aes(fill = Origin), color = "black", outlier.size = -1, lwd = 0.3) +
  geom_jitter(position=position_jitterdodge(jitter.width = 0.1), aes(color = Origin), size = 1) +
  theme_classic() +
  ylim(0,1.2) +
  theme(
    legend.position = "right",
    axis.text.y = element_text(size = 10, colour = "black"),
    axis.text.x = element_text(size = 11, colour = "black"),
    axis.title = element_text(size = 12, colour = "black", face="bold"),
    legend.title = element_text(face = "bold", size = 13),
    legend.text= element_text(size = 13)) +
  labs(x = "", y = "Frequency of NK cells") +
  scale_color_manual(values = c("#000000","#000000")) +
  scale_fill_manual(values = c("#6699FF","#FFCC33")) +
  theme(axis.text.x = element_text(angle = 45,  hjust=1)) +
  geom_pwc(
    aes(group = Origin), tip.length = 0.03,
    method = "wilcox_test",
    label.size = 4,
    label = "p.signif",
    bracket.nudge.y = -0.08, y.position = 1.2,
    symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns")))
dev.off()


# -- Alluvial Plot - which cells belong to which cluster in original tumor NK UMAP -- Fig. 2E
#load original dataset
NK <- readRDS("path/to/NK_only.RDS")

#add annotation to integrated heatly/tumor dataset
integrated$Annotation_tumor <- NA
integrated$Annotation_tumor <- NK$Fine_Annotation[match(rownames(NK@meta.data),rownames(integrated@meta.data))]
integrated$tumor <- ifelse(integrated$Annotation_tumor %in% c("NK3_CAMK4", "NK3_ENTPD1"), 
                           "tumor-associated", "non tumor-associated")

#retrieve table with 
table <- dittoBarPlot(integrated, var = "Annotation", group.by = "Annotation_tumor", retain.factor.levels = F,data.out=T)

cell_freq <- table$data
cell_freq$label <- factor(cell_freq$label, levels = rev(levels(integrated$Annotation)))
cell_freq$grouping <- factor(cell_freq$grouping, levels = levels(NK$Fine_Annotation))

cell_freq$Prevalence <- ifelse(cell_freq$label %in% c("taNK3"), "tumor-associated", "non tumor-associated")
cell_freq$Prevalence <- factor(cell_freq$Prevalence, levels = c("tumor-associated", "non tumor-associated"))

#select colors for each cluster present in both datasets
Colors <-  c(
  "NK1" = "#E15759",
  "NK2" = "#86BCB6",
  "NK3" = "#59A14F",
  "taNK3" = "#9966CC",
  "NKprol" = "#CC99FF",
  "ILC1" =  "#6699FF",
  "ILC3" =  "#CCCCCC",
  "NK1_CCL4" = "#FF9D9A",
  "NK_NFKB1" = "#499894",
  "NK3_GZMK" = "#59A14F",
  "NK3_CAMK4" = "#8CD17D",
  "NK3_ENTPD1" = "#9966CC",
  "NKprol" = "#CC99FF",
  "tumor-associated" = "#E15759",
  "non tumor-associated" = "#CCCCCC",
  "NK1_FGFBP2"= "#E15759")

cell_freq_new <- cell_freq[,c("label", "grouping", "Prevalence", "count")]
cell_freq_long <- to_lodes_form(data.frame(cell_freq_new),
                                key = "Demographic", value = "Group", id = "Cohort",
                                axes = 1:2)

cell_freq_long$Group <- factor(cell_freq_long$Group, levels = c("NK1_FGFBP2", "NK1_CCL4",
                                                                "NK1","NK2", "NK_NFKB1", "NK3",
                                                                "NK3_GZMK", "taNK3","NK3_CAMK4","NK3_ENTPD1", 
                                                                "NKprol", "ILC1", "ILC3"))


pdf("path/to/Alluvial.pdf", height = 7, width = 6) 
ggplot(data = cell_freq_long,
       aes(x = Demographic, stratum = Group, alluvium = Cohort, y = count, alpha = Prevalence)) +
  geom_alluvium(aes(fill=Prevalence, color = after_scale(fill)), width = 0.01)+
  geom_stratum(aes(fill=Group), width = 0.125, linewidth=0.8, alpha=1, color="white") +
  geom_text(stat = "stratum", aes(label = Group),
            nudge_x = -0.08,hjust = 1, data = function(x) x[x$Demographic == "label", ], size=5,  alpha = 1) +
  geom_text(stat = "stratum", aes(label = Group),
            nudge_x = +0.08,hjust = 0, data = function(x) x[x$Demographic == "grouping", ], size=5,  alpha = 1) +
  theme_classic() +
  xlab(NULL)+
  ylab("Cell count")+
  scale_fill_manual(name="Prevalence", values = Colors, breaks = c("non tumor-associated","tumor-associated")) +
  scale_alpha_manual(values = c(0.7, 1))+ 
  scale_x_discrete(labels = c("healthy and \ntumor lung", "tumor \nlung"))+
  guides(alpha="none", fill = guide_legend(override.aes = list(alpha = 1)))+
  theme(
    legend.position = "right",
    axis.line.x = element_blank(), 
    axis.line.y = element_blank(), 
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 18, colour = "black", face="bold"),
    axis.title.y = element_blank(), 
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(size = 18, colour = "black", face="bold"),
    legend.title = element_text(face = "bold", size = 13),
    legend.text= element_text(size = 13)) 
dev.off()


saveRDS("path/to/NK_HLvsTL.RDS")












