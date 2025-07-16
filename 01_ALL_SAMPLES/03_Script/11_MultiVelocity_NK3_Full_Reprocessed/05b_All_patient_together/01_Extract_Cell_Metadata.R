#In Seurat V5:
#Convert the Seurat Object to an AnnData Object
Object_To_Convert = readRDS( file.path( PATH_EXPERIMENT_RAWDATA , "RDS_FILES",SAMPLE_NK3_Only_Preprocessed) ) 

#Renaming the clusters according to new namings
levels(Object_To_Convert@meta.data[["Fine_Annotation_WNN"]]) = c("NK1_FGFBP2", "NK1_CCL4" , "NK2_NFKB1" , "NK2_XCL1" ,"NK3_GZMK" , "NK3_CAMK4" ,  "NK3_ENTPD1",  "NK_Prolif"  , "ILC1" ,  "ILC3",   "NK_DNAJB1" )


#Reset the ident
Object_To_Convert = SetIdent(Object_To_Convert, value= variable_to_sort)

#Subsetting to remove the NK_DNAJB1
DimPlot(Object_To_Convert, reduction = "WNN_UMAP")

Object_To_Convert = subset(Object_To_Convert , idents= c("NK3_GZMK" , "NK3_CAMK4" ,  "NK3_ENTPD1") )

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

if (DO_EXTRACT_EMBED==TRUE){
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

}


#save for later:
saveRDS(Object_To_Convert,file.path(PATH_ANALYSIS_OUTPUT, "subseted_seurat_object", "Subseted_seurat_Object.rds") ) 


