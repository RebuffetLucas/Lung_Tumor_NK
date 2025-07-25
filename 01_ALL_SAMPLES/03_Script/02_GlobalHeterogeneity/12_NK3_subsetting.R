library(Seurat)
library(Signac)
library(ggpubr)
library(dplyr)
set.seed(1234)

NK <- readRDS("path/to/NK_only.RDS")

# -- subset NK3 for trajectory and SCENIC+ analysis -- 
NK3 <- subset(NK, Fine_Annotation == "NK3_GZMK" |
                  Fine_Annotation == "NK3_CAMK4" |
                  Fine_Annotation == "NK3_ENTPD1")

#table(NK3$orig.ident)
#remove BS897 as there were too less cells for integration :( 
NK3 <- subset(NK3, orig.ident == "BS897", invert = T)

#RNA integration
DefaultAssay(NK3) <- "RNA"
NK3[["RNA"]] <- split(NK3[["RNA"]], f = NK3$orig.ident)
NK3 <- NormalizeData(NK3)
NK3 <- FindVariableFeatures(NK3, nfeatures = 500)

features <- VariableFeatures(NK3)

# Custom filtering using external gene blacklist
load("/scicore/home/zippeliu/serger0000/Multiome/analysis/data/exclude.gene.misc.human.v4.RData") 
f.feat <- !(features %in% all.gene.ignore.df[["seu.id"]]) &
  !(grepl("^(RP[LS]|Rp[ls])", features, perl=T)) &
  !(grepl("^MT-|^MTRN|^MTAT|^MTND|^MRP",features,perl=T)) &
  !(grepl("^RP[0-9]+-",features,perl=T)) &
  !(features %in% c('XIST')) #per hand gene rein

features.flt <-features[f.feat]

VariableFeatures(NK3) <- features.flt

NK3 <- ScaleData(NK3, features = rownames(NK3), vars.to.regress = c("percent.mt"))
NK3 <- RunPCA(NK3)

NK3 <- IntegrateLayers(object = NK3, method = CCAIntegration,
                       orig.reduction = "pca", new.reduction = "integrated.cca",
                       verbose = FALSE)


NK3[["RNA"]] <- JoinLayers(NK3[["RNA"]])

# ATAC integration - get anchors from RNA dataset
NK3_List <- SplitObject(NK3, split.by = "orig.ident")

anchors <- FindIntegrationAnchors(object.list = NK3_List,assay = c(rep("RNA", 10)),
                                  normalization.method = c("LogNormalize"),
                                  anchor.features = features.flt)

DefaultAssay(NK3) <- "peaks"
NK3 <- RunTFIDF(NK3)
NK3 <- FindTopFeatures(NK3, min.cutoff = 10)
NK3 <- RunSVD(NK3)

NK3 <- IntegrateEmbeddings(
  anchorset = anchors,
  reductions = NK3[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:25,
  verbose = T)

#WNN integration
NK3 <- FindMultiModalNeighbors(
  object =NK3_RNA,
  reduction.list = list("integrated.cca", "integrated_lsi"),
  dims.list = list(1:10, 1:15),
  verbose = TRUE)

#UMAP
NK3 <- RunUMAP(
  object = NK3,
  nn.name = "weighted.nn",
  return.model = TRUE,
  verbose = TRUE,
  reduction.name = "WNN_UMAP")

#plot previous clustering
NK3$Fine_Annotation <- droplevels(NK3$Fine_Annotation)
NK3$Fine_Annotation <- factor(NK3$Fine_Annotation, 
                              levels = c("NK3_GZMK", "NK3_CAMK4", "NK3_ENTPD1"))

Colors <- c( "#8CD17D","#FF994E","#6699FF")

pdf(file = "path/to/DimPlot_NK3.pdf",width = 3, height = 2.5)
DimPlot(NK3, reduction = "WNN_UMAP", label = T, group.by = "Fine_Annotation", cols = Colors,
        label.size =3,label.box = T, repel = T, raster = F, pt.size = 0.3) + 
  theme(plot.title = element_blank()) & NoAxes() & NoLegend()
dev.off()

























