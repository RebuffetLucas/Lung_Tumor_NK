library(Seurat)
library(Signac)
library(ggpubr)
library(dplyr)
library(UCell)
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

pdf(file = "path/to/DimPlot_NK3.pdf",width = 3, height = 2.5) #Fig. 2H
DimPlot(NK3, reduction = "WNN_UMAP", label = T, group.by = "Fine_Annotation", cols = Colors,
        label.size =3,label.box = T, repel = T, raster = F, pt.size = 0.3) + 
  theme(plot.title = element_blank()) & NoAxes() & NoLegend()
dev.off()


#enrichment with taNK sig - Fig.2H bottom-right

DEGs_taNKs <- read.csv("path/to/DEGs_taNKs.csv")

#get top30 marker
nontaNK <- rownames(taNK[order(taNK$avg_log2FC), ])[1:30]
taNK <- rownames(taNK[order(-taNK$avg_log2FC), ])[1:30]

NK3 <- AddModuleScore_UCell(NK3, features = list(taNK_sig=taNK, nontaNK_sig=nontaNK),name = "")


pdf("path/to/sigs_Vln_taNK3.pdf", width = 6, height = 3.5)
VlnPlot(NK3, c("nontaNK_sig", "taNK_sig"),group.by = "Fine_Annotation_WNN", pt.size = 0, cols = Colors) &
  ylab("UCell signature score") & 
  theme(axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 13, color = "black", face = "bold"),
        plot.title = element_blank(),
        axis.title.x = element_blank()) & NoLegend() 
dev.off()


# -- Volcano - DEGs - taNK3 vs. non taNK3 -- Fig. S4F

#UMAP colored by ta and non taNK3s
pdf("path/to/UMAP_NK3_tank.pdf", width = 3.5, height = 3.5)
DimPlot(NK3, reduction = "WNN_UMAP", label = T,
        group.by = "tumor",
        cols = c("#E15759", "#499894"),
        label.size = 4,label.box = T, repel = T, raster = F) & NoAxes() & NoLegend()
dev.off()

#Volcano
DEGs_taNK3s <- FindMarkers(NK3, ident.1 = "taNK3", ident.2 = "NK3", min.pct = 0.25,
                          test.use = "MAST", logfc.threshold = 0)

pdf("path/to/taNK3vsnontaNK3_volcano.pdf", width = 8, height = 8)
EnhancedVolcano(DEGs_taNK3s,
                lab = DEGs_taNK3s$X,
                x = 'avg_log2FC',
                y = 'p_val_adj',
                titleLabSize = 25,
                subtitleLabSize = 18,
                FCcutoff = 0.3,
                ylab = bquote(~-Log[10] ~ italic(P[adj])),
                xlim = c(-2.3, 2.3),
                pointSize = 3.0,
                legendIconSize = 3.0,
                labSize = 4,
                pCutoff = 0.0001,
                legendPosition = "bottom",
                selectLab = c("ENTPD1","GZMA", "ITGAE", "ITGA1", "RGS1", "CTLA4",
                              "CD96",  "IL32", "KLRC1", "KLRC2",
                              "PLPP1", "CSF1", "IRF4", 
                              "PRF1","GZMK", "IFNGR1", "NKG7", "CCL4", "CRTAM", "CD226", "CAMK4",
                              "PDE3B", "AHR",   "CBLB","PRF1"),
                drawConnectors = T,
                colConnectors = 'black',
                gridlines.minor = F, 
                gridlines.major = F,
                boxedLabels = TRUE,
                colAlpha = 4/5,
                col = c("#CCCCCC","#CCCCCC","#CCCCCC","#CC3333")) +
  theme(plot.title = element_blank(),
        plot.subtitle = element_blank(),
        axis.text = element_text(color = 'black'),
        axis.title = element_text(color = 'black'),
        axis.line = element_line(color = 'black'),
        legend.text = element_blank(),
        legend.position = "none")
dev.off()

# -- GSEA - taNK3 -- Fig. S4G

sigs <- readRDS( "path/to/sig.RDS")

DEGs_taNK3s <- MarkerH[order(MarkerH$avg_log2FC, decreasing = T),]

#make a named vector
gl <-DEGs_taNK3s$avg_log2FC
names(gl)<-make.names(DEGs_taNK3s$X, unique = T)

TERM2GENE <- enframe(sig, name = "term", value = "gene") %>%
  unnest(gene)

#perform GSEA
gsea_res <- GSEA(
  geneList = gl,
  TERM2GENE = TERM2GENE,
  pvalueCutoff = 0.05,
  minGSSize = 1,
  maxGSSize = 500
)

df <- as.data.frame(gsea_res)

# transform -log10(p.adjust) if you want that as color
df$log10padj <- -log10(df$p.adjust)

# Keep top N terms
df_all <- df %>% arrange(p.adjust)

# Plot
pdf("path/to/NK3_GSEA_sigs.pdf", height = 4, width = 7)
ggplot(df_all, aes(x = NES, y = reorder(Description, NES), size = setSize, color = log10padj)) +
  geom_point() +
  scale_color_viridis_c(option = "D") +
  theme_bw() +
  labs(y = "Gene Set", x = "NES",
       color = "-log10(adj. p-value)", size = "Gene Set Size") +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 70, vjust = 1, hjust = 1, size = 13, colour = "black"),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14, colour = "black"),
        axis.title.y = element_text(size = 13, colour = "black")) 
dev.off()

#saveResults
write.csv(df_all, "path/to/NK3_GSEA_dataframe.csv")


saveRDS(NK3, "path/to/NK3_only.RDS")












