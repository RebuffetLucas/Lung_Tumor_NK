#This script run the trajectory analysis for the entire subsets of Dim together
# Perform the transcriptionnal trajectories analysis (Monocle 3) on NK Cells

# Focus on Trajectory GZMK => ENTPD1
library(monocle3)
library(Seurat)
library(SeuratData)

# Version of seurat Wrapper that actually works
# remotes::install_github('satijalab/seurat-wrappers', ref = "b8feea013e7e19a46e935684b510399ffe0b6740" #
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(magrittr)
library(viridis)
library(circlize)

#Build the seurat object
#library(cerebroApp)
#Recreate relevent Seurat Object


#Params
LIST_PARAM_CONTROL =list(ncenter=107)
root_pr_nodes = "Y_7"
CLUSTER_TO_PARTITION_OUT  = "NK3_ENTPD1"

TRAJ  = "Traj_To_CAMK4" # "Traj_To_CAMK4"
RNA_OR_GENEACTIV = "REGULONS" # "GeneActivity"


NUMBER_GENES_DYNAMIC_HEATMAP_SMALL = 100 #print them all
#Load the saved slots and recreate a Seurat Object from scratch

counts = readRDS(file.path(PATH_EXPERIMENT_OUTPUT, "10_DiffusionMap", "Object_before_RUNFast_MNN",  "counts_all.rds")) # Using all the RNA

test_counts = read_csv("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Lung_Basel/01_ALL_SAMPLES/05_Output/13_SCENICplus_Analysis_NK3_Resctricted/04e_Step5_Downstream_Analysis/AUC_Tables/direct_gene_based_AUC.csv")
test_counts  = as.data.frame(test_counts)

#Reformat to match the format of counts data
#Change row names to cell names
rownames(test_counts) <- test_counts$Cell
test_counts$Cell <- NULL  # Optionally remove the Cell column now that it's rownames

# Split and reformat the row names
new_rownames <- sub(".*-(CSS\\d+)___.*", "\\1_\\1", rownames(test_counts))  # extract CSS ID
barcodes <- sub("^(.*?)-CSS\\d+___.*", "\\1", rownames(test_counts))        # extract barcode
rownames(test_counts) <- paste0(new_rownames, "-", barcodes)
rownames(test_counts) <- gsub("^([A-Za-z0-9]+)_[A-Za-z0-9]+-(.+)$", "\\1_\\2", rownames(test_counts))

#Keep only the rows present in test_counts and in counts
common_rows <- intersect(rownames(test_counts), rownames(counts))
test_counts <- test_counts[common_rows, , drop = FALSE]

# Reorder to match counts
test_counts <- test_counts[rownames(counts)[rownames(counts) %in% rownames(test_counts)], , drop = FALSE]

# Filter the regulons to keep
# Extract column names
regulon_names <- colnames(test_counts)

# Use regular expressions to extract base names and modes
pattern <- "^(.*?_direct)_(\\+\\/\\+|\\+\\/\\-|\\-\\/\\+|\\-\\/\\-)_\\(\\d+g\\)$"
regulon_df <- data.frame(
  full_name = regulon_names,
  stringsAsFactors = FALSE
)

# Extract base name and regulation type
regulon_df$base_name <- sub(pattern, "\\1", regulon_df$full_name)
regulon_df$mode <- sub(pattern, "\\2", regulon_df$full_name)

# Remove entries that didn't match (if any)
regulon_df <- regulon_df[grepl(pattern, regulon_df$full_name), ]

# Split into activators and repressors
activator_df <- regulon_df[regulon_df$mode %in% c("+/+", "+/-"), ]
repressor_df <- regulon_df[regulon_df$mode %in% c("-/+", "-/-"), ]

# Initialize final lists
filtered_activators <- c()
filtered_repressors <- c()

# Filter activators: prefer +/+ over +/-
for (tf in unique(activator_df$base_name)) {
  variants <- activator_df[activator_df$base_name == tf, ]
  if ("+/+" %in% variants$mode) {
    filtered_activators <- c(filtered_activators, variants$full_name[variants$mode == "+/+"])
  } else {
    filtered_activators <- c(filtered_activators, variants$full_name[1])
  }
}

# Filter to keep only the regulons that contain "+/+" in their name
filtered_activators_plus_plus <- filtered_activators[grepl("\\+\\/\\+", filtered_activators)]

if (Keep_only_plus_plus_regulons ==TRUE){
  filtered_activators = filtered_activators_plus_plus
}


# Filter repressors: prefer -/+ over -/-
for (tf in unique(repressor_df$base_name)) {
  variants <- repressor_df[repressor_df$base_name == tf, ]
  if ("-/+" %in% variants$mode) {
    filtered_repressors <- c(filtered_repressors, variants$full_name[variants$mode == "-/+"])
  } else {
    filtered_repressors <- c(filtered_repressors, variants$full_name[1])
  }
}


#Filter the object if necessary
if(REGULON_TO_USE == "Filtered_Activators"){
  test_counts = test_counts[,filtered_activators]
} else if(REGULON_TO_USE == "Filtered_Repressors"){
  test_counts = test_counts[,filtered_repressors]
}




#Create the object
#Load the saved slots and recreate a Seurat Object from scratch
data = readRDS( file.path(PATH_EXPERIMENT_OUTPUT, "10_DiffusionMap", "Object_before_RUNFast_MNN",  "log_norm_data_all.rds"))
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

NK_Seurat@graphs = object_for_reduction@graphs

PBMC= NK_Seurat
NK_Seurat = subset(NK_Seurat, cells = rownames(test_counts))

#Reinject the regulon data
NK_Seurat@assays[["RNA"]]@counts = t(test_counts)
NK_Seurat@assays[["RNA"]]@data = t(test_counts)

NK_Seurat = subset(NK_Seurat ,subset = Fine_Annotation_WNN_2 == CLUSTER_TO_PARTITION_OUT , invert = TRUE )
DimPlot(NK_Seurat)


cds= as.cell_data_set(NK_Seurat)
cds= estimate_size_factors(cds)


#Assign partitions
recreate.partitions <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partitions) <- cds@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)
recreate.partitions


#Assign cluster informations
list.cluster <- NK_Seurat@active.ident
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster

#Assign partition  2 to NK3

Cells_To_Partition_Out =  names(list.cluster[which(list.cluster==c(CLUSTER_TO_PARTITION_OUT))])

levels(recreate.partitions) = c(1,2)
for(name in Cells_To_Partition_Out) {
  recreate.partitions[[name]] <- factor("2", levels = c("1", "2"))
}

cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions


#Plot before trajectory
partition.before.traj <-plot_cells(cds, color_cells_by = "partition", label_groups_by_cluster = F,  group_label_size = 5) + theme(legend.position = "right")
partition.before.traj


cluster.before.traj <-plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F,  group_label_size = 5) + theme(legend.position = "right")
cluster.before.traj


#Learn Trajectory
#cds <- learn_graph(cds, use_partition = F, learn_graph_control  = list(ncenter=280), close_loop = FALSE) #Best so far
#cds <- learn_graph(cds,  use_partition = TRUE, close_loop = FALSE)
#cds <- learn_graph(cds,  use_partition = TRUE, close_loop = FALSE, learn_graph_control  = list(minimal_branch_len=30))

cds <- learn_graph(cds,  use_partition = TRUE, close_loop = FALSE, learn_graph_control=LIST_PARAM_CONTROL)


#cds <- learn_graph(cds,  use_partition = FALSE,learn_graph_control = list( nn.k=30), close_loop = FALSE)
#cds <- learn_graph(cds, learn_graph_control = list( nn.k=20), use_partition = FALSE, close_loop = FALSE)
#cds <- learn_graph(cds, learn_graph_control  = list( nn.k=15), close_loop = FALSE)
#cds <- learn_graph(cds, use_partition = F, close_loop = FALSE)

plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F ,
           label_branch_points = F, label_roots = F, label_leaves = F, cell_size = 0.5,
           graph_label_size = 5, label_principal_points = T, trajectory_graph_segment_size = 2)


plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F ,
           label_branch_points = F, label_roots = F, label_leaves = F, cell_size = 0.5,
           graph_label_size = 5, label_principal_points = F, trajectory_graph_segment_size = 2)



#Stop here to pursue with Branch analysis
#Order cells in PseudoTime
cds <- order_cells(cds, reduction_method = "UMAP", root_pr_nodes = root_pr_nodes)

plot_pseudo_time = plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = T,
           label_branch_points = F, label_roots = F, label_leaves = F)


print(plot_pseudo_time)




ggsave(
  file.path(PATH_ANALYSIS_OUTPUT , RNA_OR_GENEACTIV , "Trajectories_plots", TRAJ ,  paste0("PseudotimePlot" , TRAJ , ".pdf" )),
  plot = plot_pseudo_time,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 15,
  height = 12.5,
  units = "cm",
  dpi = 600,
  limitsize = TRUE,
  bg = NULL
)



cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))

#Visualize Pseudotime
ggplot(data.pseudo, aes(monocle3_pseudotime, Fine_Annotation_WNN_2, fill = Fine_Annotation_WNN_2 )) + geom_boxplot()
ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(Fine_Annotation_WNN_2, monocle3_pseudotime), fill = Fine_Annotation_WNN_2)) + geom_boxplot()


#Filter the most significant
#Best of the best

modulated_genes <- graph_test(cds, neighbor_graph = "principal_graph", cores = 12)

#Top50
modulated_genes %>% dplyr::arrange(q_value, desc(morans_I)) %>%
  dplyr::filter(status == "OK") %>%
  dplyr::filter(!grepl("RPS", rownames(.))) %>%
  dplyr::filter(!grepl("RPL", rownames(.))) %>%
  dplyr::filter(!grepl("MT-", rownames(.))) %>%
  dplyr::filter(q_value < 0.05 ) %>%
  head(n=NUMBER_GENES_DYNAMIC_HEATMAP_SMALL) -> top50

row.names(top50)
genes=row.names(top50)


#Order them by pseudotime
pt.pseudotime = as.data.frame(pseudotime(cds)[order(pseudotime(cds))])
#Remove cells that are not considered in this trajectory
pt.pseudotime <- pt.pseudotime[is.finite(pt.pseudotime[[1]]), , drop = FALSE]

pt.pseudotime = t(as.matrix(pt.pseudotime))

#pt.pseudotime = t(as.matrix(pseudotime(cds)[order(pseudotime(cds))]))


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
#pt.matrix <- normalized_counts(cds, norm_method = "log")[match(genes,rownames(rowData(cds))),order(pseudotime(cds))]

#Remove the cells that are not in the trajectory
#pt.matrix <- normalized_counts(cds, norm_method = "log")[match(genes,rownames(rowData(cds))),order(pseudotime(cds))]
#Keep only the cells contained n the trajectory

# Identify columns to keep (i.e., not in Cells_To_Partition_Out)
#columns_to_keep <- !(colnames(pt.matrix) %in% Cells_To_Partition_Out)

# Subset the matrix to keep only the desired columns
#pt.matrix <- pt.matrix[, columns_to_keep, drop = FALSE]
pt.matrix = as.matrix(pt.matrix)


#Create the matrix
pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
#rownames(pt.matrix) <- genes


HEIGHT = unit(9, "cm")
WIDTH = unit(6, "cm")
POLICE= gpar(fontsize = 3)


NUMBER_CLUSTS = 6


#K means with 6 groups
htkm <- Heatmap(
  pt.matrix,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),viridis(11)),
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


FILE = file.path(PATH_ANALYSIS_OUTPUT , RNA_OR_GENEACTIV , "Trajectories_plots", TRAJ , REGULON_TO_USE ,  paste0( TRAJ , "_Kmeans_Heatmap_TOP_", NUMBER_GENES_DYNAMIC_HEATMAP_SMALL, ".pdf" ))

pdf(file= FILE, width = 3, height = 20 , paper = "a4")
draw(ht_list)
dev.off()

#Ward.D2 Hierarchical Clustering
hthc <- Heatmap(
  pt.matrix,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),viridis(11)),
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


FILE = file.path(PATH_ANALYSIS_OUTPUT , RNA_OR_GENEACTIV , "Trajectories_plots", TRAJ , REGULON_TO_USE ,  paste0( TRAJ , "_WardD2_Heatmap_TOP_", NUMBER_GENES_DYNAMIC_HEATMAP_SMALL, ".pdf" ))


pdf(file= FILE, width = 3, height = 20 , paper = "a4")
draw(ht_list2)
dev.off()
