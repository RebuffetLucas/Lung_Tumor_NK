#Convert the Seurat Object to an AnnData Object
#Object_To_Convert = readRDS( file.path( PATH_EXPERIMENT_RAWDATA , "RDS_FILES",SAMPLE_NK_TUMOR) ) 


#Performing the pre-processing on patient samples one by one, lonly on cells that will be kept for MultiVelocity Analysis
# The source example is here: https://github.com/welch-lab/MultiVelo/blob/main/Examples/seurat_wnn/seurat_wnn.R

#First build the list of patients that will be used
path_velocity = file.path(PATH_EXPERIMENT_RAWDATA, "velocity", "velocity")
LIST_PATIENT =  list.dirs(path_velocity, full.names = FALSE, recursive = FALSE)

for (PATIENT in LIST_PATIENT){
#Then run the preprocessing on each of the samples:
# read in expression and accessbility data
object.data = Read10X(data.dir= file.path(path_velocity, PATIENT, paste0(PATIENT, "_filtered_feature_bc_matrix")))

#Renaming according to the new names structure
colnames(object.data[["Gene Expression"]]) = paste0(colnames(object.data[["Gene Expression"]]) , "_", PATIENT )
colnames(object.data[["Peaks"]]) = paste0(colnames(object.data[["Peaks"]]) , "_", PATIENT )


# subset for the same cells in the jointly filtered anndata object
barcodes <- read.delim(file.path(PATH_ANALYSIS_OUTPUT, "Out_01a_Filteredcells_per_indiv_sample" ,"filtered_cells",paste0(PATIENT,"_filtered_cells.txt")), header = F, stringsAsFactors = F)$V1

# preprocess RNA
seurat_object <- CreateSeuratObject(counts = object.data$`Gene Expression`[,barcodes])
seurat_object <- NormalizeData(seurat_object)
seurat_object <- FindVariableFeatures(seurat_object)
seurat_object <- ScaleData(seurat_object, do.scale = F) # not scaled for consistency with scVelo (optionally, use SCTransform)
seurat_object <- RunPCA(seurat_object, verbose = FALSE)
seurat_object <- RunUMAP(seurat_object, dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_') # optional

# preprocess ATAC
seurat_object[["ATAC"]] <- CreateAssayObject(counts = object.data$`Peaks`[,barcodes], min.cells = 1)
DefaultAssay(seurat_object) <- "ATAC"
seurat_object <- RunTFIDF(seurat_object)
seurat_object <- FindTopFeatures(seurat_object, min.cutoff = 'q0')
seurat_object <- RunSVD(seurat_object)
seurat_object <- RunUMAP(seurat_object, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_") # optional

# find weighted nearest neighbors
seurat_object <- FindMultiModalNeighbors(seurat_object, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50), k.nn = 50)
seurat_object <- RunUMAP(seurat_object, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_") # optional

DimPlot(seurat_object, reduction = "wnn.umap")

## Extract and save the neighborhood graph
file_path_save <- file.path(PATH_ANALYSIS_OUTPUT,"Out_01b_Neighborhood_Graph_Per_Indiv",PATIENT)

# Create the directory if it does not exist
if (!dir.exists(file_path_save)) {
  dir.create(file_path_save, recursive = TRUE)
  message("Directory created: ", file_path_save)
} else {
  message("Directory already exists: ", file_path_save)
}


# extract neighborhood graph
nn_idx <- seurat_object@neighbors$weighted.nn@nn.idx
nn_dist <- seurat_object@neighbors$weighted.nn@nn.dist
nn_cells <- seurat_object@neighbors$weighted.nn@cell.names

# save neighborhood graph
write.table(nn_idx,file = file.path(file_path_save,  "nn_idx.txt"), sep = ',', row.names = F, col.names = F, quote = F)
write.table(nn_dist, file = file.path(file_path_save,"nn_dist.txt"), sep = ',', row.names = F, col.names = F, quote = F)
write.table(nn_cells, file = file.path(file_path_save, "nn_cells.txt"), sep = ',', row.names = F, col.names = F, quote = F)
}

print("job finished")


if (DO_EXTRACT_EMBED == TRUE){
  
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
  
}


