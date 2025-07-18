#This script run the trajectory analysis for the entire subsets of NK3 together
# Perform the transcriptionnal trajectories analysis (Monocle 3) on NK Cells
library(monocle3)
library(Seurat)
library(SeuratData)
library(viridis)

# Version of seurat Wrapper that actually works
#remotes::install_github('satijalab/seurat-wrappers@community-vignette' ) 

library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(magrittr)

#library(cerebroApp)
#Recreate relevent Seurat Object

#Load the saved slots and recreate a Seurat Object from scratch
counts = readRDS(file.path(PATH_EXPERIMENT_OUTPUT, "10_DiffusionMap", "Object_before_RUNFast_MNN",  "counts.rds"))
data = readRDS( file.path(PATH_EXPERIMENT_OUTPUT, "10_DiffusionMap", "Object_before_RUNFast_MNN",  "log_norm_data.rds"))
metadata = readRDS(file.path(PATH_EXPERIMENT_OUTPUT, "10_DiffusionMap", "Object_before_RUNFast_MNN",  "metadata.rds"))
var.features = readRDS(file.path(PATH_EXPERIMENT_OUTPUT, "10_DiffusionMap" , "Object_before_RUNFast_MNN",  "var_features.rds"))

#Recreate a Seurat Object
NK_Seurat = CreateSeuratObject(counts = t(counts), meta.data = metadata)
NK_Seurat= NormalizeData(NK_Seurat)
NK_Seurat@assays[["RNA"]]@data  = t(data)
VariableFeatures(NK_Seurat) = var.features

#Add the Dim Reductions
object_for_reduction= readRDS(file.path(PATH_EXPERIMENT_OUTPUT, "10_DiffusionMap", "Object_before_RUNFast_MNN",  "object.rds"))

#Verify the order
sum(rownames(object_for_reduction@reductions[["WNN_UMAP"]]@cell.embeddings) == colnames(NK_Seurat))

#Inject WNN_UMAP
NK_Seurat@reductions[["UMAP"]] = object_for_reduction@reductions[["WNN_UMAP"]]
NK_Seurat  = SetIdent(NK_Seurat, value =  "Fine_Annotation_WNN_2")


#Try with the NK3 preprocessed object
#NK_Seurat  = readRDS(file.path( PATH_EXPERIMENT_RAWDATA,"RDS_FILES",SAMPLE_NK3_Only_Preprocessed ))
#NK_Seurat = readRDS(file.path(PATH_EXPERIMENT_OUTPUT, "10_DiffusionMap", "Object_before_RUNFast_MNN",  "object.rds"))

NK_Seurat@graphs = object_for_reduction@graphs

#Have a quick Look
DimPlot(NK_Seurat)
#Start from seurat object
PBMC= NK_Seurat

#Turn into cell data set
cds= as.cell_data_set(NK_Seurat)
cds= estimate_size_factors(cds)

# Not run
########
#Assign partitions
recreate.partitions <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partitions) <- cds@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)
recreate.partitions
cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions


#Assign cluster informations
list.cluster <- NK_Seurat@active.ident
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster

########
#Plot before trajectory
partition.before.traj <-plot_cells(cds, color_cells_by = "partition", label_groups_by_cluster = F,  group_label_size = 5) + theme(legend.position = "right")
partition.before.traj


cluster.before.traj <-plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F,  group_label_size = 5) + theme(legend.position = "right")
cluster.before.traj


#Learn Trajectory
#cds <- learn_graph(cds, use_partition = F, learn_graph_control  = list(ncenter=280), close_loop = FALSE) #Best so far
#cds <- learn_graph(cds,  use_partition = TRUE, close_loop = FALSE, learn_graph_control  = list(minimal_branch_len=25))

#cds <- learn_graph(cds,  use_partition = TRUE, close_loop = FALSE)


cds <- learn_graph(cds,  use_partition = TRUE, close_loop = FALSE,  learn_graph_control=LIST_PARAM_CONTROL)



#cds <- learn_graph(cds,  use_partition = TRUE, close_loop = FALSE)
#cds <- learn_graph(cds,  use_partition = FALSE,learn_graph_control = list( nn.k=30), close_loop = FALSE)
#cds <- learn_graph(cds, learn_graph_control = list( nn.k=20), use_partition = FALSE, close_loop = FALSE)
#cds <- learn_graph(cds, learn_graph_control  = list( nn.k=15), close_loop = FALSE)
#cds <- learn_graph(cds, use_partition = F, close_loop = FALSE)

plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F ,
           label_branch_points = F, label_roots = F, label_leaves = F, cell_size = 0.5,
           graph_label_size = 5, label_principal_points = T, trajectory_graph_segment_size = 2)


plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F ,
           label_branch_points = F, label_roots = F, label_leaves = F, cell_size = 0.5,
           group_label_size = 0, label_principal_points = F, trajectory_graph_segment_size = 2)

p1 = plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F ,
                label_branch_points = F, label_roots = F, label_leaves = F, cell_size = 0.5,
                group_label_size = 0, label_principal_points = F, trajectory_graph_segment_size = 2)

#Keep on with the Dynamic UMAP
#Order cells in PseudoTime
cds <- order_cells(cds, reduction_method = "UMAP", root_pr_nodes = root_pr_nodes)

plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = T,
           label_branch_points = T, label_roots = F, label_leaves = F)

p2 = plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = F ,
           label_branch_points = F, label_roots = F, label_leaves = F, cell_size = 0.5,
           graph_label_size = 5, label_principal_points = F, trajectory_graph_segment_size = 2)

p3 = plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = TRUE , label_cell_groups= TRUE,
           label_branch_points = F, label_roots = F, label_leaves = F, cell_size = 0.5, 
           group_label_size =  15, label_principal_points = F, trajectory_graph_segment_size = 1.8, trajectory_graph_color= "blue")

p4 = plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = TRUE , label_cell_groups= TRUE,
                label_branch_points = F, label_roots = F, label_leaves = F, cell_size = 0.5, 
                group_label_size =  15, label_principal_points = F, trajectory_graph_segment_size = 1.8, trajectory_graph_color= "red")



#Filter the most significant
   #Best of the best
modulated_genes <- graph_test(cds, neighbor_graph = "principal_graph", cores = 4)
genes <- row.names(subset(modulated_genes, q_value == 0 & morans_I > 0.25))
genes

  #Top50
modulated_genes %>% dplyr::arrange(q_value, desc(morans_I)) %>%
  dplyr::filter(status == "OK") %>%
  dplyr::filter(status == "OK") %>%
  dplyr::filter(q_value < 0.05 ) %>%
  head(n=50) -> top50

row.names(top50)
genes=row.names(top50)

  #Top100
modulated_genes %>% dplyr::arrange(q_value<0.05, desc(morans_I)) %>%
  dplyr::filter(status == "OK") %>%
  dplyr::filter(!grepl("RPS", rownames(.))) %>%
  dplyr::filter(!grepl("RPL", rownames(.))) %>%
  dplyr::filter(!grepl("MT-", rownames(.))) %>%
  dplyr::filter(q_value < 0.05 ) %>%
  head(n=NUMBER_GENES_DYNAMIC_HEATMAP) -> top100

row.names(top100)
genes=row.names(top100)


#Order them by pseudotime
pt.pseudotime = as.data.frame(pseudotime(cds)[order(pseudotime(cds))])
pt.pseudotime = t(pt.pseudotime)

#pt.pseudotime= as.matrix(pt.pseudotime)
#pt.pseudotime = rbind(rep(0,  length(pt.pseudotime)), pt.pseudotime)


pt.pseudotime = t(as.matrix(pseudotime(cds)[order(pseudotime(cds))]))


htpseudo <- Heatmap(
  pt.pseudotime,
  name                         = "pseudotime",
  col                          = plasma(11, alpha=1, begin = 0, end=1, direction = 1) ,
  show_row_names               = FALSE,
  show_column_names            = FALSE,
  km = 1,
  row_title_rot                = 0,
  cluster_rows                 = FALSE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
  heatmap_height =unit(0.5, "cm") )


print(htpseudo)

pt.matrix <- exprs(cds)[match(genes,rownames(rowData(cds))),order(pseudotime(cds))]
pt.matrix <- normalized_counts(cds, norm_method = "log")[match(genes,rownames(rowData(cds))),order(pseudotime(cds))]


#Can also use "normalized_counts" instead of "exprs" to use various normalization methods, for example:
#normalized_counts(cds, norm_method = "log")

pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes

HEIGHT = unit(15, "cm")
WIDTH = unit(12, "cm")
POLICE= gpar(fontsize = 3)



NUMBER_CLUSTS = 6


#K means with 6 groups
htkm <- Heatmap(
  pt.matrix,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = POLICE,
  km = NUMBER_CLUSTS,
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  show_row_dend = FALSE,
  cluster_row_slices           = TRUE,
  cluster_columns              = FALSE,
  heatmap_height =HEIGHT,
  heatmap_width = WIDTH,
  row_title = NULL)

ht_list = htpseudo %v% htkm
draw(ht_list)


FILE = file.path(PATH_ANALYSIS_OUTPUT , paste0( "TOP",NUMBER_GENES_DYNAMIC_HEATMAP, "Genes","Dynamic_Heatmap_NK3_RNA_Top100_", NUMBER_CLUSTS, "clusters",".png"))

png(file= FILE, width = 12, height = 16, units= "cm" ,  res=1200 )
draw(ht_list)
dev.off()



#Ward.D2 Hierarchical Clustering
hthc <- Heatmap(
  pt.matrix,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = POLICE,
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  show_row_dend = FALSE,
  cluster_row_slices           = TRUE,
  cluster_columns              = FALSE,
  heatmap_height =HEIGHT,
  heatmap_width = WIDTH)



ht_list2 = htpseudo %v% hthc
draw(ht_list2)


FILE = file.path(PATH_ANALYSIS_OUTPUT , paste0( "TOP",NUMBER_GENES_DYNAMIC_HEATMAP, "Genes","Dynamic_Heatmap_NK3_RNA_Top100_", "Noclusters",".png"))

png(file= FILE, width = 12, height = 16, units= "cm" ,  res=1200 )
draw(ht_list2)
dev.off()


#Save the plots of the trajectory

png(file=file.path(PATH_ANALYSIS_OUTPUT , "TrajOntheUMAP.png"), width = 30, height =25, units= "cm" ,  res=600 )
p1
dev.off()

png(file=file.path(PATH_ANALYSIS_OUTPUT , "TrajPseudotime.png"), width = 30, height =25, units= "cm" ,  res=600 )
p2
dev.off()

png(file=file.path(PATH_ANALYSIS_OUTPUT , "TrajPseudotimeblue.png"), width = 30, height =25, units= "cm" ,  res=600 )
p3
dev.off()

png(file=file.path(PATH_ANALYSIS_OUTPUT , "TrajPseudotimered.png"), width = 30, height =25, units= "cm" ,  res=600 )
p4
dev.off()




