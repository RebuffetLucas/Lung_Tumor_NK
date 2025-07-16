#Loading and preparing the data ###
#Load the saved slots and recreate a Seurat Object from scratch

if( Choose_feature_gene_activty_method == "Top_VST"){
  counts = readRDS(file.path(PATH_ANALYSIS_OUTPUT, "Object_GeneActivty_Before_RUNFast_MNN",  "count.data.gene_activity_top1000_vst.rds"))
  data = readRDS( file.path(PATH_ANALYSIS_OUTPUT, "Object_GeneActivty_Before_RUNFast_MNN",  "log.norm.data.gene_activity_top1000_vst.rds"))
  metadata = readRDS(file.path(PATH_ANALYSIS_OUTPUT, "Object_GeneActivty_Before_RUNFast_MNN",  "metadata.rds"))
  var.features = readRDS(file.path(PATH_ANALYSIS_OUTPUT, "Object_GeneActivty_Before_RUNFast_MNN",  "best_vst_genes_list_top1000_vst.rds"))
} else if( Choose_feature_gene_activty_method == "q95") {
  counts = readRDS(file.path(PATH_ANALYSIS_OUTPUT, "Object_GeneActivty_Before_RUNFast_MNN",  "count.data.gene_activity_q95.rds"))
  data = readRDS(file.path(PATH_ANALYSIS_OUTPUT, "Object_GeneActivty_Before_RUNFast_MNN",  "log.norm.data.gene_activity_q95.rds"))
  metadata = readRDS(file.path(PATH_ANALYSIS_OUTPUT, "Object_before_RUNFast_MNN",  "metadata.rds"))
  var.features = readRDS(file.path(PATH_ANALYSIS_OUTPUT, "Object_GeneActivty_Before_RUNFast_MNN",  "var.features.gene_activity_q95.rds"))
}


#Create the Seurat Object
PBMC = CreateSeuratObject(counts = t(counts), meta.data = metadata)
PBMC= NormalizeData(PBMC)
PBMC@assays[["RNA"]]@data  = t(data)
VariableFeatures(PBMC) = var.features


#Diffusion map (FastMNN)
object_d = PBMC

object_d2 <- RunFastMNN(object.list = SplitObject(object_d, 
                                                  split.by = "orig.ident"),
                        features=var.features,
                        assay = "RNA")


matrix=t(object_d2@assays$mnn.reconstructed@data) %>% as.matrix()

dm<-DiffusionMap(data=matrix,
                 censor_val = 30,
                 censor_range = c(30, 40),
                 verbose=TRUE
)

plot(dm)

#Save
saveRDS(dm , file.path(PATH_ANALYSIS_OUTPUT,"Object_GeneActivity_After_RUNFast_MNN",paste0("DiffusionMap_NK3", Choose_feature_gene_activty_method,".rds")))

dm= readRDS( file.path(PATH_ANALYSIS_OUTPUT,"Object_GeneActivity_After_RUNFast_MNN",paste0("DiffusionMap_NK3", Choose_feature_gene_activty_method,".rds")))

plot(dm)

dm_coor=dm@eigenvectors
#Try to remove DC1 and connected values

dm@eigenvectors = dm@eigenvectors[,2:20]
dm@eigenvalues = dm@eigenvalues[2:20]


#Root with automatic choice:

if (Choose_root_method == "Automatic"){
  
dpt_d <- DPT(dm)
} else if (Choose_root_method == "Supervised"){
  
#Root with NK3_GZMK as initial
group=factor(metadata[rownames(dm_coor),"Fine_Annotation_WNN_2"])
y=aggregate(dm_coor[,1],by=list(group),median)

# Sort groups by median values
y2 = setNames(y[,2], y[,1])  # Assign group names to median values
y2 = sort(y2)  # Sort in ascending order

y2=y[,2] %>% setNames(y[,1])
y2=sort(y2)

# Set the root to the highest dm2 value
root = which(dm_coor[,2] == max(dm_coor[,2]))
dpt_d <- DPT(dm, tips=root)

}

pt<- dpt_d[[paste0("DPT",root)]]
pt[pt>quantile(pt,0.99,na.rm=TRUE)]=NA
pt[pt<quantile(pt,0.01,na.rm=TRUE)]=NA

df=data.frame("pseudotime"=pt,cluster=group)

#rownames(df)=names(dm$DC1)  
rownames(df)=names(dm$DC2)  

#for colors 
pt_color<- dpt_d[[paste0("DPT",root)]]
pt_color[pt_color>quantile(pt_color,0.99,na.rm=TRUE)]=quantile(pt_color,0.99,na.rm=TRUE)
pt_color[pt_color<quantile(pt_color,0.01,na.rm=TRUE)]=quantile(pt_color,0.01,na.rm=TRUE)

clustering_d=metadata[rownames(dm_coor),"Fine_Annotation_WNN_2"] %>%
  as.character() %>% setNames(rownames(dm_coor))



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
  geom_point(shape=21, size=2) +  theme_bw() + theme(panel.grid.major = element_blank(),   panel.grid.minor = element_blank()) +
  scale_shape_identity(aes( fill = cluster, inherit.aes = TRUE))  + scale_fill_manual(values = col_firstCl)

#Version avec cercles
p2=ggplot(df2, aes( DC1, DC2, fill= pseudotime2)) +
  geom_point(aes(shape=21),size=2) +  theme_bw() + theme(panel.grid.major = element_blank(),   panel.grid.minor = element_blank()) +
  scale_shape_identity(aes(colour= "black", fill = pseudotime2, inherit.aes = TRUE))  + scale_fill_gradientn(colours = viridis::plasma(11, alpha=1, begin = 0, end=1, direction = 1))

p=p1+p2
print(p)


#Try with DC2 and DC3
p1=ggplot(df2, aes( DC2, DC3, fill= cluster)) +
  geom_point(shape=21, size=2) +  theme_bw() + theme(panel.grid.major = element_blank(),   panel.grid.minor = element_blank()) +
  scale_shape_identity(aes( fill = cluster, inherit.aes = TRUE))  + scale_fill_manual(values = col_firstCl)

#Version avec cercles
p2=ggplot(df2, aes( DC2, DC3, fill= pseudotime2)) +
  geom_point(aes(shape=21),size=2) +  theme_bw() + theme(panel.grid.major = element_blank(),   panel.grid.minor = element_blank()) +
  scale_shape_identity(aes(colour= "black", fill = pseudotime2, inherit.aes = TRUE))  + scale_fill_gradientn(colours = viridis::plasma(11, alpha=1, begin = 0, end=1, direction = 1))

p=p1+p2
print(p)

#Save them
ggsave(
  file.path(PATH_ANALYSIS_OUTPUT, "Figures", "GeneActivity", paste0("Diffusion_map_DC2_DC3","_AFTER_DC1_REMOVAL_",Choose_feature_gene_activty_method,".pdf")),
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
  file.path(PATH_ANALYSIS_OUTPUT, "Figures", "GeneActivity", paste0("Pseudotime_Plot","_AFTER_DC1_REMOVAL_" ,Choose_feature_gene_activty_method,".pdf")),
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

#Write as rds
sum(df2$pseudotime == df2$pseudotime2)
df2$pseudotime2 = NULL
saveRDS(df2, file.path(PATH_ANALYSIS_OUTPUT,  "Object_GeneActivity_After_RUNFast_MNN" ,paste0("dataframe_DC_Pseudotime_","_AFTER_DC1_REMOVAL_" ,Choose_feature_gene_activty_method,".rds")))

#write as csv for python reading
df2 = readRDS( file.path(PATH_ANALYSIS_OUTPUT,  "Object_GeneActivity_After_RUNFast_MNN" ,paste0("dataframe_DC_Pseudotime_","_AFTER_DC1_REMOVAL_" ,Choose_feature_gene_activty_method,".rds")))
df2$cellname = rownames(df2)
write_csv(df2 , file.path(PATH_ANALYSIS_OUTPUT,  "Object_GeneActivity_After_RUNFast_MNN" ,paste0("dataframe_DC_Pseudotime_","_AFTER_DC1_REMOVAL_" ,Choose_feature_gene_activty_method,".csv")))





                       