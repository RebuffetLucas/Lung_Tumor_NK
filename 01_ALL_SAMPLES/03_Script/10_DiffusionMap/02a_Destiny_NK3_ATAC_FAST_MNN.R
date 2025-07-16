#Running UMAP with FAstMNN Batchcorrection first
#Loading and preparing the data ###

# Load required libraries
library(ggplot2)
library(reshape2)
library(dplyr)
library(patchwork)

# Load objects
PBMC <- readRDS(file.path(PATH_EXPERIMENT_RAWDATA, "RDS_FILES", SAMPLE_NK3_Only_Preprocessed))

set.seed(SEED)

#Strategy 1: Extracting the TopFeatures q95
# Gene Activity Data
DefaultAssay(PBMC) <- "GeneActivity"
PBMC <- FindTopFeatures(PBMC, min.cutoff = "q95")
var.features.gene_activity_q95 <- VariableFeatures(PBMC)
count.data.gene_activity_q95 <- FetchData(object = PBMC, assay = "GeneActivity", layer = "counts", vars = var.features.gene_activity_q95)
log.norm.data.gene_activity_q95<- FetchData(object = PBMC, assay = "GeneActivity", layer = "data", vars = var.features.gene_activity_q95)


#Strategy 2: Extracting the top 1000 based on VST method
# Gene Activity Data
count.data.gene_activity <- FetchData(object = PBMC, assay = "GeneActivity", layer = "counts", vars = rownames(PBMC@assays[["GeneActivity"]]@data))
log.norm.data.gene_activity <- FetchData(object = PBMC, assay = "GeneActivity", layer = "data", vars = rownames(PBMC@assays[["GeneActivity"]]@data))

#Create a Seurat Object
PBMC2 = CreateSeuratObject(counts = PBMC@assays[["GeneActivity"]]@counts)
PBMC2 = NormalizeData(PBMC2)

#Injecting the correct slot
PBMC2@assays[["RNA"]]@layers[["data"]] =  PBMC@assays[["GeneActivity"]]@data

#Variable features top 1000
PBMC2 = FindVariableFeatures(PBMC2, nfeatures = 1000)
best_vst_genes_list = VariableFeatures(PBMC2)

#Extracting the corresponding matrixes of count and log norm data
count.data.gene_activity_top1000_vst <- FetchData(object = PBMC, assay = "GeneActivity", layer = "counts", vars = best_vst_genes_list)
log.norm.data.gene_activity_top1000_vst <- FetchData(object = PBMC, assay = "GeneActivity", layer = "data", vars = best_vst_genes_list)

#Look for basic variable features on RNA
DefaultAssay(PBMC) = "RNA"
Stored_variable_features = VariableFeatures(PBMC) #Previous variable features from full preprocessing
PBMC = FindVariableFeatures(PBMC, nfeatures = 1000)

#Have a quick look
p1 = VariableFeaturePlot(PBMC2) + ggtitle("Based on Gene Activity")
p2 = VariableFeaturePlot(PBMC) + ggtitle("Based on real RNA")

p1 + p2

#Save the figures
pdf(file.path(PATH_ANALYSIS_OUTPUT, "Verif_GeneActivity_Slot", "Standardized_Variance_Comparison_GeneActity_vs_RNA.pdf"),
    width = 24 , height = 12 )
print(p1 + p2)
dev.off()

#Compare the intersection between the two methods
INTERSECT = intersect(var.features.gene_activity_q95 , best_vst_genes_list )
print(INTERSECT)


metadata  = PBMC@meta.data

#Save the slots
saveRDS(metadata , file.path(PATH_ANALYSIS_OUTPUT, "Object_GeneActivty_Before_RUNFast_MNN",  "metadata.rds"))

#q95
saveRDS(count.data.gene_activity_q95 , file.path(PATH_ANALYSIS_OUTPUT, "Object_GeneActivty_Before_RUNFast_MNN",  "count.data.gene_activity_q95.rds"))
saveRDS(log.norm.data.gene_activity_q95 , file.path(PATH_ANALYSIS_OUTPUT, "Object_GeneActivty_Before_RUNFast_MNN",  "log.norm.data.gene_activity_q95.rds"))
saveRDS(var.features.gene_activity_q95 , file.path(PATH_ANALYSIS_OUTPUT, "Object_GeneActivty_Before_RUNFast_MNN",  "var.features.gene_activity_q95.rds"))

#vst 1000
saveRDS(count.data.gene_activity_top1000_vst , file.path(PATH_ANALYSIS_OUTPUT, "Object_GeneActivty_Before_RUNFast_MNN",  "count.data.gene_activity_top1000_vst.rds"))
saveRDS(log.norm.data.gene_activity_top1000_vst , file.path(PATH_ANALYSIS_OUTPUT, "Object_GeneActivty_Before_RUNFast_MNN",  "log.norm.data.gene_activity_top1000_vst.rds"))
saveRDS(best_vst_genes_list , file.path(PATH_ANALYSIS_OUTPUT, "Object_GeneActivty_Before_RUNFast_MNN",  "best_vst_genes_list_top1000_vst.rds"))


#Save object V1 and V2
#Diffusion map (FastMNN)
object_d = PBMC
#Remove all the non-RNA slots to make it easier
DefaultAssay(object_d) = "GeneActivity"
object_d@assays[["ATAC"]] = NULL
object_d@assays[["peaks"]] = NULL
object_d@assays[["RNA"]] = NULL
object_d@assays[["chromvar"]] = NULL

saveRDS(object_d ,  file.path(PATH_ANALYSIS_OUTPUT, "Object_GeneActivty_Before_RUNFast_MNN",  "object_GeneActivity.rds"))
saveRDS(PBMC2 ,  file.path(PATH_ANALYSIS_OUTPUT, "Object_GeneActivty_Before_RUNFast_MNN",  "object_GeneActivity_PBMC2.rds"))













#Then move to seurat_metaNK docker
splited_object = SplitObject(object_d,  split.by = "orig.ident")

object_d2 <- RunFastMNN(object.list = splited_object,
                        features=var.features, assay = "RNA")

matrix=t(object_d2@assays$mnn.reconstructed@data) %>% as.matrix()


dm<-DiffusionMap(data=matrix,
                 censor_val = 30,
                 censor_range = c(30, 40),
                 verbose=TRUE
)

plot(dm)

#Save

saveRDS(dm , file.path(PATH_ANALYSIS_OUTPUT,"DiffusionMap_NK3.rds"))

dm= readRDS( file.path(PATH_ANALYSIS_OUTPUT,"DiffusionMap_NK3.rds"))

plot(dm)
dm_coor=dm@eigenvectors

#Root with automatic choice:
dpt_d <- DPT(dm)

#Root with NK3_GZMK as initial
group=factor(metadata[rownames(dm_coor),"Fine_Annotation_WNN_2"], levels=c("NK3_GZMK","NK3_CAMK4","NK3_ENTPD1"))

y=aggregate(dm_coor[,1],by=list(group),median)
y2=y[,2] %>% setNames(y[,1])
y2=sort(y2)
if (names(y2)[1]=="NK3_GZMK") {
  root=which(dm_coor[,1]==max(dm_coor[,1]))
} else {
  root=which(dm_coor[,1]==max(dm_coor[,1]))   
}

dpt_d <- DPT(dm, tips=root)


pt<- dpt_d[[paste0("DPT",root)]]
pt[pt>quantile(pt,0.99,na.rm=TRUE)]=NA
pt[pt<quantile(pt,0.01,na.rm=TRUE)]=NA

df=data.frame("pseudotime"=pt,cluster=group)
rownames(df)=names(dm$DC1)  

#for colors 
pt_color<- dpt_d[[paste0("DPT",root)]]
pt_color[pt_color>quantile(pt_color,0.99,na.rm=TRUE)]=quantile(pt_color,0.99,na.rm=TRUE)
pt_color[pt_color<quantile(pt_color,0.01,na.rm=TRUE)]=quantile(pt_color,0.01,na.rm=TRUE)

clustering_d=metadata[rownames(dm_coor),"Fine_Annotation_WNN_2"] %>%
  as.character() %>% setNames(rownames(dm_coor))




#Other tries
#plot(dpt_d, root= 1, col= col_firstCl[    clustering_d])
#plot(dpt_d, pch=20, col= unname(col_firstCl[    clustering_d]), frame.plot = TRUE)  + theme_bw() + theme(panel.grid.major = element_blank(),   panel.grid.minor = element_blank())
#plot(dpt_d, pch=20, col= unname(col_firstCl[    clustering_d]), frame.plot = TRUE )  + theme_bw() + theme(panel.grid.major = element_blank(),   panel.grid.minor = element_blank())

#plot(dpt_d, root= 1, col= col_firstCl[    clustering_d])
#plot(dpt_d, col= col_firstCl[    clustering_d], divide = 3, dcs = c(-1,-3,2), pch = 20)
#plot(dpt_d,  divide = 3, dcs = c(-1,-3,2), pch = 20)


# 2D plots with GGPLOT
pt<- dpt_d[[paste0("DPT",root)]]

  #Reverse DC1
dm_coor[,"DC1"] = -dm_coor[,"DC1"]

object_d=subset(PBMC,cells=rownames(dm_coor))
object_d[["dm"]]=CreateDimReducObject(dm_coor, key="DC")
object_d$pt=pt

df2=data.frame(dm_coor[,1:3],pseudotime=pt,pseudotime2=pt_color)

Annotations = object_d[[Annotation_of_interest]]

df2$cluster = Annotations$Fine_Annotation_WNN_2





#######################################################
#######################################################
#Version avec cercles
p1=ggplot(df2, aes( DC1, DC2, fill= cluster)) +
  geom_point(aes(shape=21), size=2) +  theme_bw() + theme(panel.grid.major = element_blank(),   panel.grid.minor = element_blank()) +
  scale_shape_identity(aes(colour= "black", fill = cluster, inherit.aes = TRUE))  + scale_fill_manual(values = col_firstCl)

#Version avec cercles
p2=ggplot(df2, aes( DC1, DC2, fill= pseudotime2)) +
  geom_point(aes(shape=21),size=2) +  theme_bw() + theme(panel.grid.major = element_blank(),   panel.grid.minor = element_blank()) +
  scale_shape_identity(aes(colour= "black", fill = pseudotime2, inherit.aes = TRUE))  + scale_fill_gradientn(colours = viridis::plasma(11, alpha=1, begin = 0, end=1, direction = 1))

p=p1+p2
print(p)



#Save them
ggsave(
  file.path(PATH_ANALYSIS_OUTPUT, "Figures", "Diffusion_map.pdf"),
  plot = p,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 38,
  height = 15,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)


#######################################################
#######################################################


# Pseudotime Plot

p3 = ggplot(df, aes( pseudotime, cluster,col = cluster)) +  
  ggbeeswarm::geom_quasirandom( alpha=0.7,
                                cex=1,
                                show.legend = FALSE,
                                
                                groupOnX=FALSE)+ylab("")  + 
  ggtitle("Pseudotime in Dataset4") + theme_light()+
  theme_light(base_size = 25) +
  stat_summary(
    aes(group = cluster), fun = median, fun.min = median, fun.max = median,
    geom = "crossbar", color = "black", width = 0.7, lwd = 0.2,
    
    # add this bit here to your stat_summary function
    position=position_dodge(width=0.75)
  ) +
  scale_color_manual(values=col_firstCl)

p3


#Fig4e
ggsave(
  file.path(PATH_ANALYSIS_OUTPUT, "Figures", "Pseudotime_Plot.pdf"),
  plot = p3,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 30,
  height = 18,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)








                       