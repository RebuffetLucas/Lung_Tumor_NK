#Convert the Seurat Object to an AnnData Object

Object_To_Convert = readRDS( file.path( PATH_EXPERIMENT_RAWDATA , "RDS_FILES",SAMPLE_NK_TUMOR) ) 

#saving into H5ad
SaveH5Seurat(Object_To_Convert , filename = file.path( PATH_EXPERIMENT_RAWDATA , "RDS_FILES", "NK_Tumor.h5Seurat") )
Convert(file.path( PATH_EXPERIMENT_RAWDATA , "RDS_FILES", "NK_Tumor.h5Seurat"), dest = "h5ad")


#Extract and save_metadata
cell_data = Object_To_Convert@meta.data
write.csv(cell_data ,file = file.path(PATH_ANALYSIS_OUTPUT, 'cell_data.csv')  )

#remotes::install_github("mojaveazure/seurat-disk@v0.0.0.9015")
#packageVersion("SeuratDisk")


table(cell_data$sample)
DimPlot(Object_To_Convert, reduction = "WNN_UMAP")
Object_To_Convert

#Extract the embeddings:
# Create a list of matrices and their corresponding file names
tables_to_save <- list(
RNA_UMAP_TABLE = Object_To_Convert@reductions[["UMAP_RNA"]]@cell.embeddings,
ATAC_UMAP_TABLE = Object_To_Convert@reductions[["ATAC_umap"]]@cell.embeddings,
WNN_UMAP_TABLE = Object_To_Convert@reductions[["WNN_UMAP"]]@cell.embeddings
)

# Specify the file path where you want to save the CSV
file_path_save <- file.path(PATH_ANALYSIS_OUTPUT,"Embeddings")

# Loop through the list and save each matrix as a CSV file
for (name in names(tables_to_save)) {
  file_path <- file.path(file_path_save, paste0(name, ".csv"))
  print(file_path)
  # Save each matrix as a CSV file
  write.csv(tables_to_save[[name]], file = file_path, row.names = TRUE)
  # Print a message to confirm the file was saved
  cat("Saved:", file_path_save, "\n")
}


## Extract and save the neighborhood graph
file_path_save <- file.path(PATH_ANALYSIS_OUTPUT,"Neighborhood_Graph")

# extract neighborhood graph
nn_idx <- Object_To_Convert@neighbors$weighted.nn@nn.idx
nn_dist <- Object_To_Convert@neighbors$weighted.nn@nn.dist
nn_cells <- Object_To_Convert@neighbors$weighted.nn@cell.names

# save neighborhood graph
write.table(nn_idx,file = file.path(file_path_save,  "nn_idx.txt"), sep = ',', row.names = F, col.names = F, quote = F)
write.table(nn_dist, file = file.path(file_path_save,"nn_dist.txt"), sep = ',', row.names = F, col.names = F, quote = F)
write.table(nn_cells, file = file.path(file_path_save, "nn_cells.txt"), sep = ',', row.names = F, col.names = F, quote = F)

