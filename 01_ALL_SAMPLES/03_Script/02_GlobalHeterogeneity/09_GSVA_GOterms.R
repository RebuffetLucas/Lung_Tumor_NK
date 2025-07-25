library(Seurat)
library(RColorBrewer)
library(GSVA)
library(escape)
library(ComplexHeatmap)
library(circlize)

set.seed(1234)

NK <- readRDS("path/to/NK_only.RDS")

#agg again per cluster - remove ILCs

agg <- AggregateExpression(subset(NK, Fine_Annotation == "ILC1"| Fine_Annotation == "ILC3", slot = "scale.data",
                                  invert = T), assays = "RNA",return.seurat = T, group.by = c("orig.ident", "Fine_Annotation"))

agg_noseu <- AggregateExpression(subset(NK, Fine_Annotation == "ILC1"| Fine_Annotation == "ILC3", slot = "scale.data",
                                  invert = T), assays = "RNA",return.seurat = F, group.by = c("orig.ident", "Fine_Annotation"))

agg@meta.data[c("Patient", "cluster")] <- do.call(rbind, strsplit(agg$orig.ident,"_"))
agg$cluster <- gsub('-', '_', agg$cluster)

#escape package can be used to quickly load the GO terms
gene.sets_GO <- getGeneSets(library = "C5", subcategory = "GO:BP")

#now perform GSVA of GO gene sets on the pseudobulked dataset

gsva.es <- GSVA::gsva(as.matrix(agg_noseu[["RNA"]]), gene.sets_GO, verbose = TRUE)

gsva.es <- as.data.frame(gsva.es)
gsva.es <- t(gsva.es)

agg_GO <- AddMetaData(agg, gsva.es)
df <- data.frame(agg_GO[[]])

#only take the ones with FDR < 0.05
output_sig <- getSignificance(df, group = "cluster", fit = "KW")
output_sig_p <- subset(output_sig, output_sig$FDR < 0.05)

#subset dataframe by significant GO terms
df_p <- df[, rownames(output_sig_p)]

#delete non immune relevant from heatmap
non_immune_terms <- c("GOBP_CELLULAR_RESPONSE_TO_BRAIN_DERIVED_NEUROTROPHIC_FACTOR_STIMULUS","GOBP_POSITIVE_REGULATION_OF_TROPHOBLAST_CELL_MIGRATION",
                       "GOBP_POSITIVE_REGULATION_OF_SKELETAL_MUSCLE_TISSUE_GROWTH","GOBP_SENSORY_PERCEPTION_OF_SWEET_TASTE","GOBP_ENDOTHELIAL_CELL_FATE_COMMITMENT",
                       "GOBP_NOREPINEPHRINE_BIOSYNTHETIC_PROCESS","GOBP_REGULATION_OF_MESODERM_DEVELOPMENT","GOBP_WNT_SIGNALING_PATHWAY_INVOLVED_IN_HEART_DEVELOPMENT",
                       "GOBP_MESENCHYMAL_TO_EPITHELIAL_TRANSITION","GOBP_CANONICAL_WNT_SIGNALING_PATHWAY_INVOLVED_IN_HEART_DEVELOPMENT",
                       "GOBP_CANONICAL_WNT_SIGNALING_PATHWAY_INVOLVED_IN_CARDIAC_MUSCLE_CELL_FATE_COMMITMENT","GOBP_BRANCH_ELONGATION_OF_AN_EPITHELIUM",
                       "GOBP_ANTEROGRADE_DENDRITIC_TRANSPORT_OF_NEUROTRANSMITTER_RECEPTOR_COMPLEX","GOBP_REGULATION_OF_ICOSANOID_SECRETION",
                       "GOBP_REGULATION_OF_FAT_CELL_DIFFERENTIATION", "GOBP_POSITIVE_REGULATION_OF_FAT_CELL_DIFFERENTIATION",
                       "GOBP_BROWN_FAT_CELL_DIFFERENTIATION","GOBP_FAT_CELL_DIFFERENTIATION","GOBP_REGULATION_OF_MONOCYTE_EXTRAVASATION",
                       "GOBP_TYPE_I_PNEUMOCYTE_DIFFERENTIATION","GOBP_NEGATIVE_REGULATION_OF_MULTICELLULAR_ORGANISM_GROWTH",
                       "GOBP_MICROGLIAL_CELL_MEDIATED_CYTOTOXICITY","GOBP_HIPPOCAMPAL_NEURON_APOPTOTIC_PROCESS","GOBP_RESPONSE_TO_MOLECULE_OF_FUNGAL_ORIGIN",
                       "GOBP_CENTRAL_NERVOUS_SYSTEM_MATURATION","GOBP_REGULATION_OF_RETINOIC_ACID_RECEPTOR_SIGNALING_PATHWAY",
                       "GOBP_NEGATIVE_REGULATION_OF_RUFFLE_ASSEMBLY","GOBP_NEGATIVE_REGULATION_OF_MICROGLIAL_CELL_ACTIVATION",
                       "GOBP_REGULATION_OF_MICROGLIAL_CELL_ACTIVATION","GOBP_REGULATION_OF_HEPATOCYTE_APOPTOTIC_PROCESS",
                       "GOBP_NEGATIVE_REGULATION_OF_HEPATOCYTE_APOPTOTIC_PROCESS","GOBP_POSITIVE_REGULATION_OF_PROTEIN_CONTAINING_COMPLEX_DISASSEMBLY",
                       "GOBP_REGULATION_OF_RHO_PROTEIN_SIGNAL_TRANSDUCTION","GOBP_POSITIVE_REGULATION_OF_RHO_PROTEIN_SIGNAL_TRANSDUCTION",
                       "GOBP_STEROID_METABOLIC_PROCESS","OBP_KETONE_BIOSYNTHETIC_PROCESS","GOBP_MAMMARY_GLAND_EPITHELIAL_CELL_DIFFERENTIATION",
                       "GOBP_BODY_MORPHOGENESIS","GOBP_HEAD_MORPHOGENESIS","GOBP_POSITIVE_REGULATION_OF_EXTRACELLULAR_MATRIX_ASSEMBLY",
                       "GOBP_POSITIVE_REGULATION_OF_NEURON_DIFFERENTIATION","GOBP_PROTEIN_SULFATION","GOBP_REGULATION_OF_BLOOD_VESSEL_REMODELING",
                       "GOBP_ACTIVATION_OF_MEIOSIS","GOBP_REGULATION_OF_VASCULAR_WOUND_HEALING","GOBP_NEGATIVE_REGULATION_OF_VASCULAR_WOUND_HEALING",
                       "GOBP_RESPONSE_TO_AMYLOID_BETA","GOBP_REGULATION_OF_BIOMINERALIZATION","GOBP_REGULATION_OF_BONE_MINERALIZATION",
                       "GOBP_BONE_MINERALIZATION","GOBP_REGULATION_OF_OSSIFICATION","GOBP_POSITIVE_REGULATION_OF_OSSIFICATION",
                       "GOBP_NEGATIVE_REGULATION_OF_MUSCLE_TISSUE_DEVELOPMENT","GOBP_MAST_CELL_PROLIFERATION","GOBP_NEGATIVE_REGULATION_OF_RETINOIC_ACID_RECEPTOR_SIGNALING_PATHWAY",
                       "GOBP_NEGATIVE_REGULATION_OF_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_IN_RESPONSE_TO_OSMOTIC_STRESS",
                       "GOBP_PYROPTOSIS")

df_no <- df_p[,!names(df_p) %in% non_immune_terms]

mat_gsva <- t(scale(df_no))

mypalette <- rev(brewer.pal(n = 11, "RdBu"))
myBreaks <- seq(-1.45, 1.8, length.out = 11)


Colors <- c( "#E15759" , "#FF9D9A","#86BCB6","#499894","#8CD17D","#FF994E","#6699FF",
             "#CC99FF")

ht_GO = Heatmap(mat_gsva, name = "z-score",
                column_split = df$cluster,
                cluster_columns = F,
                show_column_dend = FALSE,
                cluster_column_slices = TRUE,
                column_title_gp = gpar(fontsize = 10),
                column_gap = unit(1, "mm"),
                cluster_rows =T,
                show_row_dend = FALSE,
                col=colorRamp2(myBreaks, mypalette),
                row_names_gp = gpar(fontsize = 9),
                column_title_rot = 45,
                top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = Colors))),
                show_column_names = F,
                use_raster = F,
                raster_quality = 4,
                row_names_side = "left")

pdf("path/to/HeatMap_GOterms_GSVA.pdf", width = 15, height =25)
draw(ht_GO, padding = unit(c(2, 130, 2, 2), "mm")) 
dev.off()


