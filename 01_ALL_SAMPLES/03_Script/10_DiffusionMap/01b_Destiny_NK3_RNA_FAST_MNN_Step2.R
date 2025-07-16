#Loading and preparing the data ###
#Run in Seurat_Meta_NK Docker
#Load the saved slots and recreate a Seurat Object from scratch
counts = readRDS(file.path(PATH_ANALYSIS_OUTPUT, "Object_before_RUNFast_MNN",  "counts.rds"))
data = readRDS( file.path(PATH_ANALYSIS_OUTPUT, "Object_before_RUNFast_MNN",  "log_norm_data.rds"))
metadata = readRDS(file.path(PATH_ANALYSIS_OUTPUT, "Object_before_RUNFast_MNN",  "metadata.rds"))
var.features = readRDS(file.path(PATH_ANALYSIS_OUTPUT, "Object_before_RUNFast_MNN",  "var_features.rds"))


#Create the Seurat Object
PBMC = CreateSeuratObject(counts = t(counts), meta.data = metadata)
PBMC= NormalizeData(PBMC)
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
saveRDS(dm , file.path(PATH_ANALYSIS_OUTPUT,"DiffusionMap_NK3.rds"))

dm= readRDS( file.path(PATH_ANALYSIS_OUTPUT,"DiffusionMap_NK3.rds"))

plot(dm)
dm_coor=dm@eigenvectors

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

# Set the root to the highest dm1 value
root = which(dm_coor[,1] == min(dm_coor[,1]))
dpt_d <- DPT(dm, tips=root)

}

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



# 2D plots with GGPLOT
pt<- dpt_d[[paste0("DPT",root)]]

  #Reverse DC1
#dm_coor[,"DC1"] = -dm_coor[,"DC1"]

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


#Write as rds
sum(df2$pseudotime == df2$pseudotime2)
df2$pseudotime2 = NULL
saveRDS(df2, file.path(PATH_ANALYSIS_OUTPUT,  "dataframe_DC_Pseudotime.rds"))


#write as csv for python reading
df2 = readRDS( file.path(PATH_ANALYSIS_OUTPUT,  "dataframe_DC_Pseudotime.rds"))
df2$cellname = rownames(df2)
write_csv(df2 , file.path(PATH_ANALYSIS_OUTPUT,  "dataframe_DC_Pseudotime.csv"))







                       