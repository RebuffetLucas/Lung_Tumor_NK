library(chromVAR)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Seurat)
library(Signac)
library(dplyr)
library(viridis)
library(ComplexHeatmap)
library(circlize)
set.seed(123456)

NK <- readRDS("path/to/NK_only.RDS")
DefaultAssay(NK) <- "peaks"

#runChromvar

NK <- RunChromVAR(
  object = NK,
  genome = BSgenome.Hsapiens.UCSC.hg38)


#find diff. accessible motifs for fine clusters

NK <- SetIdent(NK, value = "Fine_Annotation")

Motifs <- FindAllMarkers(object= NK, logfc.threshold = 0,
                          assay = "chromvar", mean.fxn = rowMeans, 
                          fc.name = "avg_diff")

write.csv(Motifs, "to/path/Chromvar_TFs.csv")

#aggregate chromvar assay + plot Heatmap

Ch <- AggregateExpression(NK, 
                          assay = "chromvar", 
                          return.seurat = F, 
                          group.by = c("Fine_Annotation"))

JASPAR24 <- readRDS("path/to/motifsTF_DF_JASPAR2024.RDS") #dataframe of motif name with respective TF name

#convert Motifs into TF names
chr_df <- as.data.frame(Ch$chromvar)
chr_df$TF <- JASPAR24$TF[match(JASPAR24$Motif, rownames(chr_df))]
chr_df_u <- chr_df[!duplicated(chr_df$TF), ]
rownames(chr_df_u) <- chr_df_u$TF
chr_df_u$TF <- NULL

mat <- chr_df_u %>% as.matrix()

#Plotting - Fig. S2F
mypalette <- viridis(16)
myBreaks <- seq(-1.5,1.5, length.out = 16)

#selected TF as Marker Motifs per cluster
to_plot <- c("TBX2", "TBX3", "EOMES", "TBX21", "TBX20", "TBX18", "TBR1", "TBX15", 
             "ELK3", "ETV5", "ETV3", "NFYA", "CTCF", "TFEC", "BHLHE41", "BHLHE40", 
             "NFKB1", "REL",
             "TCF7","MYC",
             "NR4A1", "NR4A2", "NR2F2",
             "FOSL2::JUN", "FOSL1::JUNB", "FOS::JUND", "FOS", "FOSL1", "FOS::JUN","FOS::JUNB",
             "FOSL2::JUND", "FOSL2::JUNB","BACH1", "JDP2", "BNC2", "BATF", "BATF::JUN", "BATF3",
             "SNAI3",
             "ZEB1", "FIGLA", "TCF4", "TCF12", "TCF3", "SNAI1", "MYOD1", "BCL11A", 
             "IRF8", "IRF7", "IRF9", "IRF4","NR2C1",
             "RORA", "RORC")

mat <- scale(t(mat[to_plot, ]))
rownames(mat) <- c("NK1_FGFBP2","NK1_CCL4","NK_NFKB1","NK2","NK3_GZMK",
                    "NK3_CAMK4","NK3_ENTPD1","NKprol" ,"ILC1","ILC3") 

pdf("path/to/Heatmaps_Chr_columnZ.pdf", width = 12, height = 3)
Heatmap(mat, name = "column scaled\nMotif activity",
        cluster_rows = F,
        show_column_dend = FALSE,
        cluster_row_slices = FALSE,
        row_title_gp = gpar(fontsize = 12),
        cluster_columns = F,
        show_row_dend = FALSE,
        col = colorRamp2(myBreaks, mypalette),
        column_names_gp = gpar(fontsize = 10),
        column_names_side = "bottom",
        column_names_rot = 45,
        show_column_names = T,
        show_row_names = T,
        row_names_side = "left",
        use_raster = F,
        border = T,
        border_gp = gpar(col = "black", lwd = 0.5))
dev.off()

saveRDS(NK, "path/to/NK_only.RDS")















