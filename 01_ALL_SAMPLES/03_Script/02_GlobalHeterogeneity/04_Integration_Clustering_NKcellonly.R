library(Seurat)
library(Signac)
library(ggpubr)
library(dittoSeq)
set.seed(1234)

NK_only <- readRDS("path/to/NK_only.RDS")

# -- RNA integration --

#Seurat v5

DefaultAssay(NK_only) <- "RNA"
NK_only[["RNA"]] <- split(NK_only[["RNA"]], f = NK_only$orig.ident)
NK_only <- NormalizeData(NK_only)
NK_only <- FindVariableFeatures(NK_only, nfeatures = 1000)

features <- VariableFeatures(NK_only)

# Custom filtering using external gene blacklist
load("path/to/exclude.gene.misc.human.v4.RData")
f.feat <- !(features %in% all.gene.ignore.df[["seu.id"]]) &
  !(grepl("^(RP[LS]|Rp[ls])", features, perl=T)) &
  !(grepl("^MT-|^MTRN|^MTAT|^MTND|^MRP",features,perl=T)) &
  !(grepl("^RP[0-9]+-",features,perl=T)) &
  !(features %in% c('XIST')) 

features.flt <-features[f.feat]
VariableFeatures(NK_only) <- features.flt

# Scale data with regression
NK_only <- ScaleData(NK_only, features = rownames(NK_only), vars.to.regress = c("percent.mt"))

#PCA
NK_only <- RunPCA(NK_only)

# Integrate using CCA
NK_only <- IntegrateLayers(object = NK_only, method = CCAIntegration,
                            orig.reduction = "pca", new.reduction = "integrated.cca",
                            verbose = FALSE)

# Collapse layers back into a single assay
NK_only[["RNA"]] <- JoinLayers(NK_only[["RNA"]])

# -- ATAC integration --

#get anchors using RNA dataset

NK_only_List <- SplitObject(NK_only, split.by = "orig.ident")
anchors <- FindIntegrationAnchors(object.list = NK_only_List, assay = "RNA",
                                  normalization.method = "LogNormalize", 
                                  anchor.features = features.flt)

# ATAC Integration

DefaultAssay(NK_only) <- "peaks"
NK_only <- RunTFIDF(NK_only)
NK_only <- FindTopFeatures(NK_only, min.cutoff = 10)
NK_only <- RunSVD(NK_only)

NK_only <- IntegrateEmbeddings(
  anchorset = anchors,
  reductions = NK_only[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:25,
  verbose = T)

# -- WNN Integration -- 

NK_only <- FindMultiModalNeighbors(
  object = NK_only,
  reduction.list = list("integrated.cca", "integrated_lsi"),
  dims.list = list(1:15, 1:20),
  verbose = TRUE)

#UMAP
NK_only <- RunUMAP(
  object = NK_only,
  nn.name = "weighted.nn",
  return.model = TRUE,
  verbose = TRUE,
  reduction.name = "WNN_UMAP")

NK_only <- FindClusters(NK_only, graph.name = "wsnn", verbose = T, 
                   resolution = seq(0, 1, 0.1))

#plot different resolutions
plt <- list()
y <- seq(0.1,1,by=0.1)

for (i in 1:10){
  plt[[i]] <- DimPlot(NK_only, reduction = "WNN_UMAP", label = T, label.size = 6,
                      group.by = paste("wsnn_res.", y[[i]], sep = "")) + ggtitle(paste("resolution", y[[i]], sep = " ")) &
    NoLegend() & NoAxes()
}

pdf(file = "path/to/res_wnn.pdf",
    height = 9, width = 12)
ggarrange(plotlist = plt, nrow = 3,ncol = 4)
dev.off()

#Coarse Clusters 0.2 - split by Patient - Fig. S1E
Colors_WNN_coarse <- c("#E15759","#FFCC00","#86BCB6","#59A14F","#999999", "#CCCCCC")

pdf(file = "path/to/DimPlot_splitPat.pdf",width = 12, height = 4)
DimPlot(NK_only, reduction = "WNN_UMAP", label = F, group.by = "wsnn_res.0.2", split.by = "orig.ident",
        cols = Colors_WNN_coarse,
        label.size = 3,label.box = T, repel = T, raster = F, ncol = 6) & NoAxes() 
dev.off()

# split by Patient - Fig. S1F
pdf(file = "path/to/Barplot_patients.pdf",width = 7, height = 5)
dittoBarPlot(NK_only, color.panel = Colors_WNN_coarse, "wsnn_res.0.2", group.by = "orig.ident",
             ylab = "Frequency of cells") &
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 13, color = "black"),
        axis.title.y=element_blank(), 
        plot.title = element_blank(), 
        legend.text = element_text(size = 12)) & coord_flip() 
dev.off()

# -- Stress cluster exclusion --

# plot stress signature - Fig. S1G -  stress score as defined in https://doi.org/10.1016/j.cell.2023.07.034 

stress_score <- c("BAG3","CALU", "DNAJB1","DUSP1", "EGR1", "FOS", "FOSB", 
                 "HIF1A", "HSP90AA1", "HSP90AB1", "HSP90B1", "HSPA1A", "HSPA1B", "HSPA6", "HSPB1", 
                 "HSPH1", "IER2", "JUN", "JUNB", "NFKBIA", "NFKBIZ", "RGS2","SLC2A3", 
                 "SOCS3", "UBC", "ZFAND2A", "ZFP36","ZFP36L1")

NK_only <- AddModuleScore_UCell(NK_only, features = list(stress_score=stress_score), name = "")

pdf(file = "path/to/VlnPlot_stressScore.pdf",width = 3, height = 2)
VlnPlot(NK_only, "stress_score", group.by =  "wsnn_res.0.2", pt.size = 0, cols = Colors_WNN_coarse) +
  ylab("Stress score") +
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 13, color = "black", face = "bold"),
        plot.title = element_blank(),
        axis.title.x = element_blank()) & NoLegend() 
dev.off()

pdf(file = "path/to/FeaturePlot_stressScore.pdf", width = 6, height = 4)
FeaturePlot(NK_only, c("HSPA1B", "HSPA1A","DNAJB1","BAG3", "HSPD1", "HSPH1"),
            cols = c("lightgrey", "#3333CC"), reduction = "WNN_UMAP",
            ncol = 3, min.cutoff = 0) & NoAxes() & NoLegend()
dev.off()

##subset stressed cluster 5
NK_only <- subset(NK_only, 
                  wsnn_res.0.2 == 1, 
                 invert = T)

#WNN UMAP
Colors_WNN <- c("#E15759","#86BCB6","#59A14F","#999999", "#CCCCCC")

p3 <- DimPlot(NK_only, reduction = "WNN_UMAP", label = F, group.by = "wsnn_res.0.2",cols = Colors_WNN,
              label.size = 3,label.box = T, repel = T, raster = F) & NoAxes() & NoLegend()

# -- FINAL RNA, ATAC, and WNN UMAPs -- 

#RNA UMAP

DefaultAssay(NK_only) <- "RNA"
ElbowPlot(NK_only)
NK_only <- FindNeighbors(NK_only, dims = 1:15, reduction = "integrated.cca")
NK_only <- RunUMAP(NK_only, dims = 1:15, reduction.name = "UMAP_RNA", reduction = "integrated.cca")
NK_only <- FindClusters(NK_only ,resolution = seq(0, 1, 0.1)) #check different resolutions

#plot different resolutions
plt <- list()
y <- seq(0.1,1,by=0.1)

for (i in 1:10){
  plt[[i]] <- DimPlot(NK_only, reduction = "UMAP_RNA", label = T, label.size = 6,
                      group.by = paste("RNA_snn_res.", y[[i]], sep = "")) + ggtitle(paste("resolution", y[[i]], sep = " ")) &
    NoLegend() & NoAxes()
}

pdf(file = "path/to/res_RNA.pdf",
    height = 9, width = 12)
ggarrange(plotlist = plt, nrow = 3,ncol = 4)
dev.off()

#resolution 0.2
NK_only$Coarse_Clustering_RNA <- NK_only$RNA_snn_res.0.2

Colors_RNA <- c("#E15759","#59A14F","#86BCB6","#999999", "#CCCCCC", "#CC99FF")

p2 <- DimPlot(NK_only, reduction = "UMAP_RNA", label = F, group.by = "Coarse_Clustering_RNA", cols = Colors_RNA,
              label.size = 3,label.box = T, repel = F, raster = F) & NoAxes() & NoLegend()


#ATAC UMAP

DefaultAssay(NK_only) <- "peaks"
DepthCor(NK_only, reduction = "integrated_lsi") #dims look all ok

NK_only <- FindNeighbors(object = NK_only, reduction = 'integrated_lsi', dims = 1:20)
NK_only <- FindClusters(object = NK_only, verbose = FALSE, graph.name = "peaks_snn", algorithm = 3,
                   resolution = seq(0, 1, 0.1))
NK_only <- RunUMAP(object = NK_only, reduction = 'integrated_lsi', dims =1:20, reduction.name = "ATAC_umap")

plt <- list()
y <- seq(0.1,1,by=0.1)

for (i in 1:10){
  plt[[i]] <- DimPlot(NK_only, reduction = "ATAC_umap", label = T, label.size = 6,
                      group.by = paste("peaks_snn_res.", y[[i]], sep = "")) + ggtitle(paste("resolution", y[[i]], sep = " ")) &
    NoLegend() & NoAxes()
}

pdf(file = "path/to/res_ATAC.pdf",
    height = 9, width = 12)
ggarrange(plotlist = plt, nrow = 3,ncol = 4)
dev.off()

#resolution 0.3

Colors_ATAC <- c("#E15759","#86BCB6","#59A14F", "#8CD17D", "#6699FF","#999999","#CCCCCC")

p1 <- DimPlot(NK_only, reduction = "ATAC_umap", label = F, group.by = "Annotation_Peaks_RNA_2", cols = Colors_ATAC,
              label.size = 3,label.box = T, repel = T, raster = F) & NoAxes() & NoLegend()


#plotting Fig. 1A
pdf(file = "path/to/DimPlots_combined.pdf",width = 8, height = 3)
p1+p2+p3 & NoLegend()
dev.off()


saveRDS(NK_only, "path/to/NK_only.RDS")












