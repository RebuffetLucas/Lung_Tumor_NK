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
library(dplyr)
library(tidyr)
library(ggalluvial)
set.seed(1234)

NK <- readRDS("path/to/NK_only.RDS")

#Fig. S3A
Colors_RNA <- c( "#E15759" ,"#FF9D9A","#499894","#8CD17D","#6699FF",
                 "#999999", "#CCCCCC")

pdf(file =  "path/to/UMAP_RNA_AnnoRNA.pdf",height = 4, width = 4)
DimPlot(NK, reduction = "UMAP_RNA", label = T,
        group.by = "Annotation_RNA",cols = Colors_RNA,
        label.size = 4,label.box = T, repel = F, raster = F) & NoAxes() & NoLegend()
dev.off()

#Fig. S3B
Colors_ATAC <- c( "#E15759" , "#FF9D9A","#499894","#8CD17D","#FF994E","#6699FF",
                  "#999999", "#CCCCCC")

pdf(file = "path/to/UMAP_ATAC_numbers.pdf",height = 4, width = 5)
DimPlot(NK, reduction = "ATAC_umap", label = T,
        group.by = "peaks_snn_res.0.5",
        cols = Colors_ATAC,
        label.size = 4,label.box = T, repel = F, raster = F) & NoAxes() & NoLegend()
dev.off()

#Fig. S3C
Colors_WNN <- c( "#E15759" , "#FF9D9A","#86BCB6","#499894","#8CD17D","#FF994E","#6699FF",
                 "#CC99FF",  "#999999", "#CCCCCC")


pdf(file = "path/to/UMAP_WNN_AnnoWNN.pdf",height = 4, width = 5)
DimPlot(NK, reduction = "WNN_UMAP", label = T,
        group.by = "Fine_Annotation_WNN",cols = Colors_WNN,
        label.size = 4,label.box = T, repel = T, raster = F) & NoAxes() & NoLegend()
dev.off()

#Alluvial Plot
anno_df<- NK@meta.data[,c("Fine_Annotation_WNN", "Annotation_RNA", "peaks_snn_res.0.5")]

df_long <- anno_df %>%
  mutate(event = row_number()) %>%
  pivot_longer(-event, names_to = "axis", values_to = "group") %>%
  count(axis, group, event) %>%
  group_by(axis) %>%
  mutate(Fraction = n / sum(n)) %>%
  ungroup()

df_long$axis <- factor(df_long$axis, levels = c("Annotation_RNA", "peaks_snn_res.0.5", "Fine_Annotation_WNN"), labels = c("RNA clustering", "ATAC clustering", "Integrated WNN clustering"))

df_long$group <- factor(df_long$group, levels = c("1","3","5","0","2","6","4","7",
                                                  "NK1_FGFBP2","NK1_CCL4",   "NK1_CX3CR1", "NK_NFKB1","NK2", "NK3_GZMK", "NK3_CAMK4",  "NK3_ENTPD1", "NKprol", "ILC1","ILC3"))

#Fig. S3D
pdf("path/to/plots/Annotation_Alluvial.pdf", height = 6, width = 6)
ggplot(df_long,
       aes(x = axis, stratum = group, alluvium = event,
           y = Fraction, fill = group)) +
  geom_flow() +
  geom_stratum(width = 0.5) +
  scale_fill_manual(values = c(Colors_RNA, Colors_WNN, Colors_ATAC)) +
  theme_minimal()+
  geom_text(
    stat = "stratum",
    aes(label = group),
    size = 3,
    color = "black"
  ) +
  scale_x_discrete(expand = c(0, 0.5))+ NoLegend() + xlab(NULL)
dev.off()

# -- CORRELATION HEATMAP --

# -- RNA-ATAC CORRELATION --

#Calculate top markers per cluster
#RNA
NK <- SetIdent(NK, value = NK$Annotation_RNA)
Markers_RNA <- FindAllMarkers(NK, assay = "RNA", min.pct = 0.25, only.pos = T) # test.use = "MAST",

#write.csv(Markers_RNA, "Marker_RNA.csv")

Markers_RNA %>%
  group_by(cluster) %>%
  arrange(cluster, p_val_adj) %>%
  slice_head(n = 50) %>%
  ungroup() -> top50_RNA

top50_RNA$cluster <- factor(top50_RNA$cluster, levels = levels(NK$Annotation_RNA))

NK_aggro_RNA <- AggregateExpression(NK, assays = "RNA", return.seurat = T, group.by = "Annotation_RNA")

NK_aggro_RNA$Annotation_RNA <- gsub("-", "_", NK_aggro_RNA$Annotation_RNA)
NK_aggro_RNA$Annotation_RNA <- factor(NK_aggro_RNA$Annotation_RNA, levels = levels(NK$Annotation_RNA))


#ATAC
NK <- SetIdent(NK, value = NK$peaks_snn_res.0.5)
Markers_ATAC <- FindAllMarkers(NK, assay = "GeneActivity", min.pct = 0.25, only.pos = T) # test.use = "MAST",

#write.csv(Markers_RNA, "Marker_RNA.csv")

Markers_ATAC %>%
  group_by(cluster) %>%
  arrange(cluster, p_val_adj) %>%
  slice_head(n = 50) %>%
  ungroup() -> top50_ATAC

top50_ATAC$cluster <- factor(top50_ATAC$cluster, levels = levels(NK$peaks_snn_res.0.5))

NK_aggro_ATAC <- AggregateExpression(NK, assays = "GeneActivity", return.seurat = T, group.by = "peaks_snn_res.0.5")

NK_aggro_ATAC$peaks_snn_res.0.5 <- gsub("g", "", NK_aggro_ATAC$peaks_snn_res.0.5)
NK_aggro_ATAC$peaks_snn_res.0.5 <- factor(NK_aggro_ATAC$peaks_snn_res.0.5, levels = levels(NK$peaks_snn_res.0.5))


gene_list <- split(top50_RNA$gene, top50_RNA$cluster)

NK_aggro_ATAC <- AddModuleScore_UCell(NK_aggro_ATAC, features = gene_list, name = "")

mat<- NK_aggro_ATAC@meta.data[,names(gene_list)] %>% as.matrix()

rownames(mat) <- gsub("g", "", rownames(mat))

mat<- scale(t(scale(mat)))

quantile(mat, c(0.05, 0.95))
mat[mat>1.8] <- 1.8
mat[mat< -1.8] <- -1.8

# Fig. S3E
pdf("path/to/plots/Corr_Heatmap_RNA_ATAC.pdf", height = 4, width = 4.4)
Heatmap(mat, name = "Gene activity\n(z-score)",
        column_split = NK_aggro_ATAC@meta.data$peaks_snn_res.0.5,
        row_split = NK_aggro_RNA@meta.data$Annotation_RNA,
        cluster_columns = F,
        cluster_column_slices = F,
        cluster_rows = F,
        #column_title_gp = gpar(fontsize = 11),
        column_gap = unit(0, "mm"),
        row_gap = unit(0, "mm"),
        show_row_dend = FALSE,
        col=rev(brewer.pal(n = 9, name = "RdBu")),
        row_names_gp = gpar(fontsize = 11),
        column_title_rot = 0,
        column_names_rot = 0,
        column_names_centered = T,
        column_names_side = "bottom",
        column_title = NULL,
        row_title = NULL,
        right_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = c(Colors_RNA)))),
        bottom_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = c(Colors_ATAC)))),
        show_column_names = T,
        show_row_names = T,
        use_raster = TRUE,
        raster_quality = 4,
        row_names_side = "right",
        rect_gp = gpar(col = "black", lwd = 3)
)
dev.off()

# -- ATAC WNN CORRELATION --

gene_list <- split(top50_ATAC$gene, top50_ATAC$cluster)

NK_aggro_WNN <- AggregateExpression(NK, assays = "GeneActivity", return.seurat = T, group.by = "Fine_Annotation_WNN")
NK_aggro_WNN$Fine_Annotation_WNN <- gsub("-", "_", NK_aggro_WNN$Fine_Annotation_WNN)
NK_aggro_WNN$Fine_Annotation_WNN <- factor(NK_aggro_WNN$Fine_Annotation_WNN, levels = levels(NK$Fine_Annotation_WNN))

NK_aggro_WNN <- AddModuleScore_UCell(NK_aggro_WNN, features = gene_list, name = "", assay = "GeneActivity")

mat<- NK_aggro_WNN@meta.data[,names(gene_list)] %>% as.matrix()
rownames(mat) <- gsub("-", "_", rownames(mat))

mat<- scale(t(scale(mat)))

quantile(mat, c(0.05, 0.95))
mat[mat>1.8] <- 1.8
mat[mat< -1.8] <- -1.8

#Fig. S3F
pdf("path/to/plots/Corr_Heatmap_ATAC_WNN.pdf", height = 4.5, width = 5) 
ht <- Heatmap(mat, name = "Gene activity\n(z-score)",
              column_split = NK_aggro_WNN@meta.data$Fine_Annotation_WNN,
              row_split = NK_aggro_ATAC@meta.data$peaks_snn_res.0.5,
              cluster_columns = F,
              cluster_column_slices = F,
              cluster_rows = F,
              column_gap = unit(0, "mm"),
              row_gap = unit(0, "mm"),
              show_row_dend = FALSE,
              col=rev(brewer.pal(n = 9, name = "RdBu")),
              row_names_gp = gpar(fontsize = 11),
              column_title_rot = 0,
              column_names_rot = 45,
              #column_names_centered = F,
              column_names_side = "bottom",
              column_title = NULL,
              row_title = NULL,
              right_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = c(Colors_ATAC)))),
              bottom_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = c(Colors_WNN)))),
              show_column_names = T,
              show_row_names = T,
              use_raster = TRUE,
              raster_quality = 4,
              row_names_side = "right",
              rect_gp = gpar(col = "black", lwd = 3)
)
draw(ht , padding = unit(c(2, 20, 2, 2), "mm"))
dev.off()

# -- RNA WNN CORRELATION --

gene_list <- split(top50_RNA$gene, top50_RNA$cluster)

NK_aggro_WNN <- AggregateExpression(NK, assays = "RNA", return.seurat = T, group.by = "Fine_Annotation_WNN")
NK_aggro_WNN$Fine_Annotation_WNN <- gsub("-", "_", NK_aggro_WNN$Fine_Annotation_WNN)
NK_aggro_WNN$Fine_Annotation_WNN <- factor(NK_aggro_WNN$Fine_Annotation_WNN, levels = levels(NK$Fine_Annotation_WNN))

NK_aggro_WNN <- AddModuleScore_UCell(NK_aggro_WNN, features = gene_list, name = "", assay = "RNA")

mat<- NK_aggro_WNN@meta.data[,names(gene_list)] %>% as.matrix()
rownames(mat) <- gsub("-", "_", rownames(mat))

#rownames(mat) <- gsub("g", "", rownames(mat))

mat<- scale(t(scale(mat)))

quantile(mat, c(0.05, 0.95))
mat[mat>1.8] <- 1.8
mat[mat< -1.8] <- -1.8


#Fig. S3G
pdf("path/to/plots/Corr_Heatmap_RNA_WNN.pdf", height = 4.5, width = 5.4)
ht <- Heatmap(mat, name = "z-score",
              column_split = NK_aggro_WNN@meta.data$Fine_Annotation_WNN,
              row_split = NK_aggro_RNA@meta.data$Annotation_RNA,
              cluster_columns = F,
              cluster_column_slices = F,
              cluster_rows = F,
              column_gap = unit(0, "mm"),
              row_gap = unit(0, "mm"),
              show_row_dend = FALSE,
              col=rev(brewer.pal(n = 9, name = "RdBu")),
              row_names_gp = gpar(fontsize = 11),
              column_title_rot = 0,
              column_names_rot = 45,
              #column_names_centered = F,
              column_names_side = "bottom",
              column_title = NULL,
              row_title = NULL,
              right_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = c(Colors_RNA)))),
              bottom_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = c(Colors_WNN)))),
              show_column_names = T,
              show_row_names = T,
              use_raster = TRUE,
              raster_quality = 4,
              row_names_side = "right",
              rect_gp = gpar(col = "black", lwd = 3)
)
draw(ht , padding = unit(c(2, 20, 2, 2), "mm"))
dev.off()

saveRDS(NK, "path/to/NK_only.RDS")








