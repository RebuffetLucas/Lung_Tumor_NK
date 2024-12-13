## @knitr DiagnosticAfterLabelTransfer

cat("# Diagnostic After Label Transfer {.tabset .tabset-fade} \n\n")

cat("## Diag heterogeneity quality Data {.tabset .tabset-fade} \n\n")

cat("### Data and Meta NK data raw \n\n")


  #Global UMAP
p01 = DimPlot(PBMC_Meta_NK_V2, cols = palette2) 

p02 = DimPlot(PBMC_Blood_All, group.by = QUERY_CLUSTER_COLNAME) 

p03 = DimPlot(PBMC_Blood_All, group.by =QUERY_ORIG_IDENT_COLNAME) 



cat("\n\n")
print(p01)
cat("\n\n")

cat("\n\n")
print( p02 + p03)
cat("\n\n")


cat("### Diag data Counts and Features \n\n")
  #Diag of Counts an Features

PBMC_Blood_All$nCount = colSums(x = PBMC_Blood_All, slot = "counts")  # nCount_RNA
PBMC_Blood_All$nFeature = colSums(x = GetAssayData(object = PBMC_Blood_All, slot = "counts") > 0)  # nFeatureRNA

p04 = VlnPlot(PBMC_Blood_All, features= "nCount", group.by = QUERY_CLUSTER_COLNAME) & stat_summary(fun.data=data_summary,color="black")
p05 = VlnPlot(PBMC_Blood_All, features= "nCount", group.by = QUERY_ORIG_IDENT_COLNAME) & stat_summary(fun.data=data_summary,color="black")

#p06 = VlnPlot(PBMC_Blood_All, features= "nCount", group.by = "tissue_Firstclust")

cat("\n\n")
print(p04)
cat("\n\n")

cat("\n\n")
print(p05)
cat("\n\n")

#cat("\n\n")
#print(p06)
#cat("\n\n")



cat("## Vizualization of the reference {.tabset .tabset-fade} \n\n")

  #Have a look at the reference built with the 4 Dataset of Meta-NK

p1 <- DimPlot(data_integrated, reduction = "umap", group.by = "Dataset")
p2 <- DimPlot(data_integrated, reduction = "umap", group.by = "FirstClust", label = TRUE, 
              repel = TRUE, cols = palette2) + NoLegend()

cat("\n\n")
print( p1 + p2 )

cat("\n\n")

p4 <- DimPlot(data_integrated, reduction = "umap", split.by = "Dataset", group.by = "FirstClust", label = TRUE, 
              repel = TRUE, cols = palette2) + NoLegend()

print(p4)
cat("\n\n")


cat("## Assessing quality of Label Transfer {.tabset .tabset-fade} \n\n")

cat("### Look at the common markers for the different populations {.tabset .tabset-fade} \n\n")

data.query = SetIdent(data.query, value = "predicted.id")

cat("#### Main pops based on V2 \n\n")

  #Main pops
Markers_NK_123 = readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/VF1/DEG/Top20_Different_Levels/TOP20_MainClusters.rds")
  
  #ViolinPlot
#Remove the very little identified clusters that could mislead dataviz in the scaled plots
if (DO_REMOVE_VERY_LITTLE_IDENT_SUBSET==TRUE){
  # Get the counts of each predicted.id
  predicted_counts <- table(data.query$predicted.id)
  # Retrieve the names with counts less than TRESHOLD_NEW_ID_REMOVING
  predicted_less_than_treshold <- names(predicted_counts[predicted_counts < TRESHOLD_NEW_ID_REMOVING])
  data.query2 = subset(data.query, subset= predicted.id %in% predicted_less_than_treshold , invert=TRUE)
}else{data.query2 = data.query
}


for (i in names(Markers_NK_123)){
  data.query2 = AddModuleScore(data.query2, features = as.list(Markers_NK_123[[i]]), pool= NULL ,name= i , seed=19)
}

p_Vln1 = VlnPlot(data.query2, features = paste0(names(Markers_NK_123),"1") , group.by = "predicted.id", cols= palette2, pt.size = 0)  & stat_summary(fun.data=data_summary,color="black")

cat("\n\n")
print(p_Vln1)
cat("\n\n")


  #Direct Dot Plot
p5 = DotPlot(data.query2, features = Markers_NK_123$NK1, cols = "RdBu") + theme(axis.text.x = element_text(angle = 90)) + ggtitle("NK1 genes") + theme(axis.text = element_text(size = 20)) 
p6 = DotPlot(data.query2, features = Markers_NK_123$NK2, cols = "RdBu") + theme(axis.text.x = element_text(angle = 90)) + ggtitle("NK2 genes") + theme(axis.text = element_text(size = 20)) 
p7 = DotPlot(data.query2, features = Markers_NK_123$NK3, cols = "RdBu") + theme(axis.text.x = element_text(angle = 90)) + ggtitle("NK3 genes") + theme(axis.text = element_text(size = 20)) 

cat("\n\n")
print(p5)
cat("\n\n")
print(p6)
cat("\n\n")
print(p7)
cat("\n\n")

cat("#### Main pops based on CITEseq \n\n")

#Main pops
Markers_NK_123_CITEseq = readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/12_CITEseq_Analysis/DEG/All_Markers_CITEseq3clusters.rds")

Markers_NK_123_CITEseq %>%
  filter(avg_log2FC>0) %>%
  filter( p_val_adj < 5e-2) %>%
  filter(!grepl("RPL|RPS|MT-", gene)) %>% #We use RPS and RPL for scoring but nor for plotind DotPlots
  group_by(cluster) %>%
  arrange(-avg_log2FC) %>%
  top_n(n = NUMBER_TOP_SCORING, wt = avg_log2FC) -> top_All

top_All %>% filter(cluster== "NK_1") -> Top_NK1
top_All %>% filter(cluster== "NK_2") -> Top_NK2
top_All %>% filter(cluster== "NK_3") -> Top_NK3


#ViolinPlot
#Remove the very little identified clusters that could mislead dataviz in the scaled plots
if (DO_REMOVE_VERY_LITTLE_IDENT_SUBSET==TRUE){
  # Get the counts of each predicted.id
  predicted_counts <- table(data.query$predicted.id)
  # Retrieve the names with counts less than TRESHOLD_NEW_ID_REMOVING
  predicted_less_than_treshold <- names(predicted_counts[predicted_counts < TRESHOLD_NEW_ID_REMOVING])
  data.query2 = subset(data.query, subset= predicted.id %in% predicted_less_than_treshold , invert=TRUE)
}else{data.query2 = data.query
}

data.query2 = AddModuleScore(data.query2, features = list(Top_NK1$gene), pool= NULL ,name= "NK1" , seed=19)
data.query2 = AddModuleScore(data.query2, features = list(Top_NK2$gene), pool= NULL ,name= "NK2" , seed=19)
data.query2 = AddModuleScore(data.query2, features = list(Top_NK3$gene), pool= NULL ,name= "NK3" , seed=19)


p_Vln2 = VlnPlot(data.query2, features = c("NK11", "NK21", "NK31") , group.by = "predicted.id", cols= palette2, pt.size = 0)  & stat_summary(fun.data=data_summary,color="black")

cat("\n\n")
print(p_Vln2)
cat("\n\n")

p8 = DotPlot(data.query2, features = Top_NK1$gene, cols = "RdBu") + theme(axis.text.x = element_text(angle = 90)) + ggtitle("NK1 genes") + theme(axis.text = element_text(size = 20)) 
p9 = DotPlot(data.query2, features = Top_NK2$gene, cols = "RdBu") + theme(axis.text.x = element_text(angle = 90)) + ggtitle("NK2 genes") + theme(axis.text = element_text(size = 20)) 
p10 = DotPlot(data.query2, features = Top_NK3$gene, cols = "RdBu") + theme(axis.text.x = element_text(angle = 90)) + ggtitle("NK3 genes") + theme(axis.text = element_text(size = 20)) 

cat("\n\n")
print(p8)
cat("\n\n")
print(p9)
cat("\n\n")
print(p10)
cat("\n\n")


cat("#### Secondary pops \n\n")
  #Secondary pops
Markers_NK_6pops = readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/VF1/DEG/Top20_Different_Levels/TOP11_6_Clusters.rds")

#Remove the very little identified clusters that could mislead dataviz in the scaled plots
if (DO_REMOVE_VERY_LITTLE_IDENT_SUBSET==TRUE){
  # Get the counts of each predicted.id
  predicted_counts <- table(data.query$predicted.id)
  # Retrieve the names with counts less than TRESHOLD_NEW_ID_REMOVING
  predicted_less_than_treshold <- names(predicted_counts[predicted_counts < TRESHOLD_NEW_ID_REMOVING])
  data.query2 = subset(data.query, subset= predicted.id %in% predicted_less_than_treshold , invert=TRUE)
}else{data.query2 = data.query
}

for (i in names(Markers_NK_6pops)){
  data.query2 = AddModuleScore(data.query2, features = as.list(Markers_NK_6pops[[i]]), pool= NULL ,name= i , seed=19)
}

p_Vln3 = VlnPlot(data.query2, features = paste0(names(Markers_NK_6pops),"1") , group.by = "predicted.id", cols= palette2, pt.size = 0)  & stat_summary(fun.data=data_summary,color="black") & stat_summary(fun.data=data_summary,color="black")

cat("\n\n")
print(p_Vln3)
cat("\n\n")

p11 = DotPlot(data.query2, features = Markers_NK_6pops$NK1A , cols = "RdBu") + theme(axis.text.x = element_text(angle = 90)) + ggtitle("NK1A genes") + theme(axis.text = element_text(size = 20)) 
p12 = DotPlot(data.query2, features = Markers_NK_6pops$NK1B, cols = "RdBu") + theme(axis.text.x = element_text(angle = 90)) + ggtitle("NK1B genes") + theme(axis.text = element_text(size = 20)) 
p13 = DotPlot(data.query2, features = Markers_NK_6pops$NK1C, cols = "RdBu") + theme(axis.text.x = element_text(angle = 90)) + ggtitle("NK1C genes") + theme(axis.text = element_text(size = 20)) 
p14 = DotPlot(data.query2, features = Markers_NK_6pops$NKint, cols = "RdBu") + theme(axis.text.x = element_text(angle = 90)) + ggtitle("NKint genes") + theme(axis.text = element_text(size = 20)) 
p15 = DotPlot(data.query2, features = Markers_NK_6pops$NK2, cols = "RdBu") + theme(axis.text.x = element_text(angle = 90)) + ggtitle("NK2 genes") + theme(axis.text = element_text(size = 20)) 
p16 = DotPlot(data.query2, features = Markers_NK_6pops$NK3, cols = "RdBu") + theme(axis.text.x = element_text(angle = 90)) + ggtitle("NK3 genes") + theme(axis.text = element_text(size = 20)) 

cat("\n\n")
print(p11)
cat("\n\n")
print(p12)
cat("\n\n")
print(p13)
cat("\n\n")
print(p14)
cat("\n\n")
print(p15)
cat("\n\n")
print(p16)
cat("\n\n")



cat("### Reliability of the prediction across predicted id \n\n")

p17 = VlnPlot(data.query2 , features= "prediction.score.max", group.by = "predicted.id" , cols = palette2, pt.size = 0) & stat_summary(fun.data=data_summary,color="black")

data.query2$predicted.id = as.factor(data.query2$predicted.id)




for(pop_predicted in levels(data.query2$predicted.id)){
  data.query2_subset = subset(data.query2, subset = predicted.id == pop_predicted)
  Vlnp_loop= VlnPlot(data.query2_subset , features= "prediction.score.max", group.by = QUERY_CLUSTER_COLNAME , pt.size = 0) & stat_summary(fun.data=data_summary,color="black") & ggtitle(pop_predicted)
  Ridge_loop= RidgePlot(data.query2_subset , features= "prediction.score.max", group.by = QUERY_CLUSTER_COLNAME )  & ggtitle(pop_predicted)
  
  cat("\n\n")
  print(Vlnp_loop)
  cat("\n\n")
  
  cat("\n\n")
  print(Ridge_loop)
  cat("\n\n")
  
  Vlnp_loop2= VlnPlot(data.query2_subset , features= "prediction.score.max", group.by = QUERY_ORIG_IDENT_COLNAME , pt.size = 0) & stat_summary(fun.data=data_summary,color="black") & ggtitle(pop_predicted)
  Ridge_loop2= RidgePlot(data.query2_subset , features= "prediction.score.max", group.by = QUERY_ORIG_IDENT_COLNAME )  & ggtitle(pop_predicted)
  
  cat("\n\n")
  print(Vlnp_loop2)
  cat("\n\n")
  
  cat("\n\n")
  print(Ridge_loop2)
  cat("\n\n")

}


cat("\n\n")
print(p17)
cat("\n\n")

p17b = VlnPlot(data.query2 , features= "prediction.score.max", group.by = QUERY_CLUSTER_COLNAME, pt.size = 0)  

cat("\n\n")
print(p17b)
cat("\n\n")

p17c = VlnPlot(data.query2 , features= "prediction.score.max", group.by = "predicted.id", split.by = QUERY_CLUSTER_COLNAME)  

cat("\n\n")
print(p17c)
cat("\n\n")

p17d = VlnPlot(data.query2 , features= "prediction.score.max", group.by = QUERY_CLUSTER_COLNAME, split.by = "predicted.id" , cols = c(palette2[1],palette2[3] ,palette2[4],palette2[2]))  
p17d

cat("\n\n")
print(p17d)
cat("\n\n")





cat("### Look at the spontaneous markers for the different populations predicted \n\n")
cat("#####" ,"DotPlot spontaneous markers","\n")
if (DO_REMOVE_VERY_LITTLE_IDENT_SUBSET==TRUE){
  # Get the counts of each predicted.id
  predicted_counts <- table(data.query$predicted.id)
  # Retrieve the names with counts less than TRESHOLD_NEW_ID_REMOVING
  predicted_less_than_treshold <- names(predicted_counts[predicted_counts < TRESHOLD_NEW_ID_REMOVING])
  data.query2 = subset(data.query, subset= predicted.id %in% predicted_less_than_treshold , invert=TRUE)
}else{data.query2 = data.query
}

data.query2 = SetIdent(data.query2, value = "predicted.id")

All_Markers = FindAllMarkers(data.query2 , only.pos = FINDMARKERS_ONLYPOS, method= FINDMARKERS_METHOD , min.pct =  FINDMARKERS_MINPCT, logfc.threshold = FINDMARKERS_LOGFC_THR , verbose = TRUE)

#Look at markers
All_Markers %>%
  group_by(cluster) %>%
  filter(!grepl("RPL|RPS|MT-", gene)) %>%
  top_n(n = FINDMARKERS_SHOWTOP, wt = avg_log2FC) -> top10


p18 = DotPlot(data.query2, features =  unique(top10$gene) , cols = "RdBu") + theme(axis.text.x = element_text(angle = 90)) + theme(axis.text = element_text(size = 20)) 

cat("\n\n")
print(p18)
cat("\n\n")

#interactive markers table
#Make the interactive table
DEG_sample <- All_Markers
DEG_sample_filtered =(DEG_sample[DEG_sample$p_val_adj<FDRCUTOFF,])
DEG_sample_filtered <- DEG_sample_filtered[order(DEG_sample_filtered$avg_log2FC ,decreasing = T),]

cat("#####" ,"Markers Table","\n")

Nb_markers=length(unique(data.query2@active.ident))
mypalette= hue_pal()(Nb_markers)
mypalette=palette2 #Vérifier si ça marche ou si ça fait bugger
print( htmltools::tagList(DT::datatable(DEG_sample_filtered,rownames = FALSE,extensions = 'Buttons', 
                                        options = list(dom = 'Blfrtip', 
                                                       buttons = c('excel', "csv"), fixedHeader = TRUE)
)%>% 
  DT::formatStyle(
    'cluster',
    backgroundColor = DT::styleEqual(sort(unique(data.query2@active.ident)),  mypalette[1:Nb_markers])
  )))


for(clust in sort(unique(Idents(data.query2)))){
  cat("\n\n")
  cat("##### Cluster ", clust,  "{.tabset .tabset-fade}\n\n")
  cat("\n\n")
  DEG_clust <- DEG_sample_filtered[DEG_sample_filtered$cluster %in% clust,]
  
  for(gene in head(DEG_clust$gene)){
    cat("\n\n")
    cat("######", gene)
    cat("\n\n")
    print(VlnPlot(data.query2, group.by = "predicted.id", features = gene, pt.size = 0))
    cat("\n\n")
  }
}



cat("## Look at the results of label transfer {.tabset .tabset-fade} \n\n")

cat("\n\n")
p22 = ggplot(data.query2@meta.data, aes(x=SecondClust, fill= predicted.id)) + geom_bar(position="fill") + theme(text = element_text(size=20),axis.text.x = element_text(angle = 90))  + scale_fill_manual(values= palette2) + ggtitle("Data Only") + theme(text = element_text(size = 20))

cat("\n\n")
#print(p22)

cat("\n\n")
p21 = ggplot(data.query2@meta.data, aes_string(x=QUERY_CLUSTER_COLNAME, fill= "predicted.id")) + geom_bar(position="fill") + theme(text = element_text(size=20),axis.text.x = element_text(angle = 90))  + scale_fill_manual(values= palette2) + ggtitle("Data Only") + theme(text = element_text(size = 20))

cat("\n\n")
print(p21)


cat("\n\n")
p23 = ggplot(PBMC_Meta_NK_V2@meta.data, aes(x=orig.ident, fill= FirstClust)) + geom_bar(position="fill") + theme(text = element_text(size=20),axis.text.x = element_text(angle = 90)) + scale_fill_manual(values= palette2) + ggtitle("Meta NK per Sample") + theme(text = element_text(size = 20))
print(p23)
cat("\n\n")

cat("\n\n")
p24 = ggplot(PBMC_Meta_NK_V2@meta.data, aes(x=Dataset, fill= FirstClust)) + geom_bar(position="fill") + theme(text = element_text(size=20),axis.text.x = element_text(angle = 90)) + scale_fill_manual(values= palette2)  + ggtitle("Meta NK per Dataset") + theme(text = element_text(size = 20))
print(p24)
cat("\n\n")

cat("## Look at the results of label transfer in the initial UMAP {.tabset .tabset-fade} \n\n")
#Reinject the prdicted.id in the initial object
predicted_id_to_add = readRDS(file.path(PATH_ANALYSIS_OUTPUT, "predicted_id_4pops.rds"))
names(predicted_id_to_add) <- sub("^Query_", "", names(predicted_id_to_add))

# Ensure the row names match
common_names <- intersect(rownames(PBMC_Blood_All@meta.data), names(predicted_id_to_add))

# Add the "predicted.id" column
PBMC_Blood_All@meta.data$predicted.id <- NA  # Initialize the column with NA
PBMC_Blood_All@meta.data[common_names, "predicted.id"] <- predicted_id_to_add[common_names]
p31 = DimPlot(PBMC_Blood_All, group.by = "predicted.id", pt.size = 0.6)

cat("\n\n")
print(p31)
cat("\n\n")


#p27 = ggplot(PBMC_Blood_All@meta.data, aes_string(x="RNA_snn_res.0.2", fill= "predicted.id")) + geom_bar(position="fill") + theme(text = element_text(size=20),axis.text.x = element_text(angle = 90))  + scale_fill_manual(values= palette2) + ggtitle("Data Only") + theme_minimal() + theme(text = element_text(size = 20))


p21 = ggplot(PBMC_Blood_All@meta.data, aes_string(x=QUERY_CLUSTER_COLNAME, fill= "predicted.id")) + geom_bar(position="fill") + theme(text = element_text(size=20),axis.text.x = element_text(angle = 90))  + scale_fill_manual(values= palette2) + ggtitle("Data Only") + theme_minimal() + theme(text = element_text(size = 20))

p21 <- ggplot(PBMC_Blood_All@meta.data, 
              aes_string(x = QUERY_CLUSTER_COLNAME, fill = "predicted.id")) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = palette2) +
  ggtitle("Data Only") +
  theme_minimal() +  # Use a minimal theme for a cleaner look
  theme(
    text = element_text(size = 20),       # General text size
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),  # Rotate x-axis labels
    plot.title = element_text(hjust = 0.5, face = "bold"),  # Center-align and bold title
    legend.title = element_text(face = "bold"),  # Bold legend title
    legend.position = "right"  # Place legend on the right
  )


cat("\n\n")
print(p21)
cat("\n\n")

#Save the figures


pdf(file.path(PATH_ANALYSIS_OUTPUT_FIGURES, "UMAP_Predicted_Id.pdf"),  width = 10, height = 6)
p31
dev.off()

pdf(file.path(PATH_ANALYSIS_OUTPUT_FIGURES, "Barplot_Predicted_Id.pdf"),  width = 10, height = 6)
p21
dev.off()

DimPlot(PBMC_Blood_All , group.by = "RNA_snn_res.0.2")

#pdf(file.path(PATH_ANALYSIS_OUTPUT_FIGURES, "Barplot_oldclust_Predicted_Id.pdf"),  width = 10, height = 6)
#p27
#dev.off()

saveRDS(PBMC_Blood_All , RDS_DIR_OUTPUT)


