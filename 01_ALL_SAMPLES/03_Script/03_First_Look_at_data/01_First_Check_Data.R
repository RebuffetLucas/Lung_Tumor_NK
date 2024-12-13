#This codes aims at having a quick Look at the data

proj <- readRDS( file.path(PATH_EXPERIMENT_RAWDATA , "RDS_FILES", SAMPLE_NK_TUMOR ))

DimPlot(proj, group.by = "Minor_Annotation")
DimPlot(proj, group.by = "Major_Annotation")


#Have a look at the data structure
proj@assays[["ATAC"]]@counts[1:10, 1:10]
proj@assays[["peaks"]]@counts[1:10, 1:10]
proj@assays[["GeneActivity"]]@counts[1:10, 1:10]
proj@assays[["chromvar"]][1:10, 1:10]

dim(proj@assays[["ATAC"]]@counts)
dim(proj@assays[["peaks"]]@counts)



intersect(rownames(proj@assays[["ATAC"]]@counts) , rownames(proj@assays[["peaks"]]@counts))

#Some tests

num_zero_rows <- sum(rowSums(proj@assays[["ATAC"]]@counts != 0) == 0)

# Print the result
num_zero_rows

help("CreateChromatinAssay")


rownames(proj@assays[["ATAC"]]@counts)[1:10]
rownames(proj@assays[["peaks"]]@counts)[1:10]
