#Try extracting the kNN graph only on the cells retained after QC
seurat_object = readRDS(file.path(PATH_EXPERIMENT_RAWDATA , "RDS_FILES", SAMPLE_NK3_Only_Preprocessed))

#List of cells to extract
cells_names =  read.delim(file.path(PATH_ANALYSIS_OUTPUT , "Out_01a_Filteredcells_per_indiv_sample" , "filtered_cells", "filtered_cells.txt"), header = FALSE)


# Reformat the names
formatted_names <- gsub("(.+)-1_(CSS\\d+)", "\\2_\\1-1", cells_names$V1)


#Focus on cells of interest
seurat_object = subset(seurat_object, cells = formatted_names )

#Re calculate the neighbors graph after subsetting
seurat_object = FindMultiModalNeighbors(seurat_object, reduction.list =  list("integrated.cca", "integrated_lsi"),
                                        dims.list = list(1:15, 1:20),
                                        verbose = TRUE)

file_path_save <- file.path(PATH_ANALYSIS_OUTPUT,"Neighborhood_Graph")

#Extract and save the neigbor graph:
# extract neighborhood graph
nn_idx <- seurat_object@neighbors$weighted.nn@nn.idx
nn_dist <- seurat_object@neighbors$weighted.nn@nn.dist
nn_cells <- seurat_object@neighbors$weighted.nn@cell.names

# save neighborhood graph
write.table(nn_idx,file = file.path(file_path_save,  "nn_idx.txt"), sep = ',', row.names = F, col.names = F, quote = F)
write.table(nn_dist, file = file.path(file_path_save,"nn_dist.txt"), sep = ',', row.names = F, col.names = F, quote = F)
write.table(nn_cells, file = file.path(file_path_save, "nn_cells.txt"), sep = ',', row.names = F, col.names = F, quote = F)

