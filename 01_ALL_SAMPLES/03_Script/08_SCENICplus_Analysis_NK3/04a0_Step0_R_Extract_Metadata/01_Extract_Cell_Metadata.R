#In Seurat V5:
#Convert the Seurat Object to an AnnData Object
Object_To_Convert = readRDS( file.path( PATH_EXPERIMENT_RAWDATA , "RDS_FILES",SAMPLE_NK_TUMOR) ) 

#saving into H5ad has to be run only once
#SaveH5Seurat(Object_To_Convert , filename = file.path( PATH_EXPERIMENT_RAWDATA , "RDS_FILES", "NK_Tumor.h5Seurat") )
#Convert(file.path( PATH_EXPERIMENT_RAWDATA , "RDS_FILES", "NK_Tumor.h5Seurat"), dest = "h5ad")

#Renaming the clusters according to new namings
levels(Object_To_Convert@meta.data[["Fine_Annotation_WNN"]]) = c("NK1_FGFBP2", "NK1_CCL4" , "NK2_NFKB1" , "NK2_XCL1" ,"NK3_GZMK" , "NK3_CAMK4" ,  "NK3_ENTPD1",  "NK_Prolif"  , "ILC1" ,  "ILC3",   "NK_DNAJB1" )

#Reset the ident
Object_To_Convert = SetIdent(Object_To_Convert, value= "Fine_Annotation_WNN")

#Subsetting to remove the NK_DNAJB1
DimPlot(Object_To_Convert, reduction = "WNN_UMAP")

Object_To_Convert = subset(Object_To_Convert , idents= "NK_DNAJB1", invert = TRUE )

DimPlot(Object_To_Convert, reduction = "WNN_UMAP")

  
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

#save for later:
saveRDS(Object_To_Convert,file.path(PATH_ANALYSIS_OUTPUT, "subseted_seurat_object", "Subseted_seurat_Object.rds") ) 

#Reload
Object_To_Convert= readRDS( file.path(PATH_ANALYSIS_OUTPUT, "subseted_seurat_object", "Subseted_seurat_Object.rds"))


#Find Neigbhors as they are lacking
Object_To_Convert = FindNeighbors(Object_To_Convert)
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





