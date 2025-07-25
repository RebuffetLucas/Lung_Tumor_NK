library(Seurat)
library(Signac)
library(UCell)
library(ggpubr)
set.seed(1234)

merged <- readRDS("path/to/merged_subset.RDS")

#---Split Object and perform SCTransform for Normalization---

NK_list <- SplitObject(merged, split.by = "orig.ident")

NK_list <- lapply(
  NK_list,
  function(obj) {
    DefaultAssay(obj) <- "RNA"
    SCTransform(
      object = obj,
      assay = "RNA",
      verbose = TRUE,
      vst.flavor = "v2",
      vars.to.regress = c("percent.mt", "percent.rb"),
      return.only.var.genes = FALSE)
  }
)

# ---Integrate Object over Patients---

#--RNA--
# Set default assay for all samples to SCT
NK_list <- lapply(NK_list, function(obj) {
  DefaultAssay(obj) <- "SCT"
  obj
})

# Select integration features
features <- SelectIntegrationFeatures(NK_list)

# Load list of genes to exclude (external file) -  exclude genes irrelevant for clustering
load("path/to/exclude.gene.misc.human.v4.RData")

# Filter features: exclude ignored genes and ribosomal genes
features.flt <- features[
  !(features %in% all.gene.ignore.df[["seu.id"]]) &
    !grepl("^(RP[LS]|Rp[ls])", features, perl = TRUE)]

# Prep for SCT-based integration
NK_list <- PrepSCTIntegration(object.list = NK_list, anchor.features = features.flt)

# Find integration anchors
anchors <- FindIntegrationAnchors(
  object.list = NK_list,
  normalization.method = "SCT",
  anchor.features = features.flt
)

# Integrate RNA dataset
integrated_NK <- IntegrateData(anchors, normalization.method = "SCT")

#--ATAC--
DefaultAssay(integrated_NK) <- "peaks"
integrated_NK <- RunTFIDF(integrated_NK)
integrated_NK <- FindTopFeatures(integrated_NK, min.cutoff = 10)
integrated_NK <- RunSVD(integrated_NK)

# integrate LSI embeddings
integrated_NK <- IntegrateEmbeddings(
  anchorset = anchors,
  reductions = integrated_NK[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:30,
  verbose = T)

# ---WNN integration ---

#PCA
DefaultAssay(integrated_NK) <- "integrated"
integrated_NK  <- RunPCA(integrated_NK, npcs = 30, verbose = T)

#WNN
integrated_NK <- FindMultiModalNeighbors(
  object = integrated_NK,
  reduction.list = list("pca", "integrated_lsi"),
  dims.list = list(1:20, 2:30),
  verbose = TRUE
)

#UMAP to remove CD3

integrated_NK <- RunUMAP(
  object = integrated_NK,
  nn.name = "weighted.nn",
  return.model = TRUE,
  verbose = TRUE,
  reduction.name = "wnn_umap")

#Clustering

integrated_NK <- FindClusters(integrated_NK, graph.name = "wsnn", algorithm = 3, verbose = T,
                              resolution = 0.2)

#plot clusters - Fig. S1B
pdf(file = "path/to/UMAP.pdf", width = 3, height = 3)
DimPlot(integrated_NK, label = T) & NoAxes() & NoLegend()
dev.off()

#calculate Markers
library(tidyverse)
integrated_NK <- PrepSCTFindMarkers(object = integrated_NK)


Markers <- FindAllMarkers(integrated_NK, min.pct = 0.25, logfc.threshold = 0.25,
                          assay = "SCT", test.use = "MAST")

Markers  %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) -> top20

#Identification of Bcell, Monotype, Tregs (cluster 9, 13, 12 and 2 based on marker genes)

#removal of T cells
T_noCD3E <- c("CD3D", "CD3G")
T_noCD3E_l <- list(T_noCD3E=T_noCD3E)

integrated_NK <- AddModuleScore_UCell(integrated_NK,
                                     assay = "SCT",
                                     features = T_noCD3E_l)

#plotting of T cell signature score - Fig. S1C

pdf("path/to/hist_CD3sig_beforesubset.pdf", width = 2, height = 1.5)
ggplot(integrated_NK@meta.data, aes(x = T_noCD3E_UCell)) +
  geom_histogram(binwidth = 0.05, fill = "cornflowerblue", color = "black", size = 0.1) +
  geom_vline(xintercept = 0.05, color = "red", linetype = "dashed", size = 0.5) +
  labs(x = "T cell signature score \n (CD3D, CD3G)", y = "Cell count") +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 7, face = "bold"),
    axis.title.y = element_text(size = 7, face = "bold"),
    axis.text.x = element_text(size = 4, color = "black"),
    axis.text.y = element_text(size = 4, color = "black"),
    panel.grid = element_blank()
  )
dev.off()

#plotting of Dotplot - Fig. S1D

pdf(file = "path/to/DotPlot.pdf", width =7, height = 3.5)
DotPlot(integrated_NK,
        features =  c("NKG7", "FGFBP2", "FCGR3A","GZMB","GLYN","KLRF1","NCAM1", "XCL1", "AREG", "C1orf21", 
                      "CD3E", "CD3G", "CD3D", "CD8A", "TRGC2","CD4","IL2RA", "CTLA4", "FOXP3",  "RORA", "CD6", "IL7R",
                      "AHR", "RORC", "IL1R1","CD79A" ,"MZB1", "JCHAIN", "CST3", "FCER1A", "XCR1", "LYZ", "CD86"),
        assay= "SCT", cols = c("RdBu")) + RotatedAxis() +
  theme(axis.text.y = element_text(size = 9),
        axis.text.x = element_text(size = 9, face = "bold"),
        legend.title = element_text(size=9),
        legend.text = element_text(size=8),
        legend.key.size = unit(0.3, "cm"),
        axis.title = element_blank()) 
dev.off()

#subset to only obtain NK cells:
# - T cells (based on CD3 signature)
# - Unwanted clusters (identified as B cells, Mono, Tregs)

NK_only <- subset(
      integrated_NK,
      subset = T_noCD3E_UCell < 0.05 &
              wsnn_res.0.2 != 9 &
              wsnn_res.0.2 != 13 &
              wsnn_res.0.2 != 12 &
              wsnn_res.0.2 != 2)

saveRDS(NK_only, "path/to/NK_only.RDS")





















