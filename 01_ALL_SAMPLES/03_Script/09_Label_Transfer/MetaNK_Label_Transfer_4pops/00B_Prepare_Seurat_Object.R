## Needs to be run only once with Seurat v5 docker.
## This code prepares a Seurat Object compatible and ready for Label transfer

# Load Query
PBMC_Blood_All = readRDS(file.path(PATH_EXPERIMENT_RAWDATA, "RDS_FILES",SAMPLE_NK_TUMOR))

#Remove the unused slots
PBMC_Blood_All@assays[["ATAC"]] = NULL
PBMC_Blood_All@assays[["peaks"]] = NULL
PBMC_Blood_All@assays[["GeneActivity"]] = NULL
PBMC_Blood_All@assays[["chromvar"]] = NULL

gc()

DefaultAssay(PBMC_Blood_All) = "RNA"
#Updating the names
levels( PBMC_Blood_All@meta.data[["Fine_Annotation_WNN"]]) = c("NK1_FGFBP2" ,"NK1_CCL4"  , "NK_NFKB1" ,  "NK2"  , "NK3_GZMK",    "NK3_CAMK4" ,  "NK3_ENTPD1" , "NK_MKI67",   "ILC1"      , "ILC3" ,   "NK_DNAJB1" )

count.data <- GetAssayData(object = PBMC_Blood_All[["RNA"]], layer = "counts")

New_object = CreateSeuratObject(counts = count.data, meta.data = PBMC_Blood_All@meta.data)

New_object@reductions[["umap"]]= PBMC_Blood_All@reductions[["WNN_UMAP"]]

New_object = SetIdent(New_object, value = "Fine_Annotation_WNN")
DimPlot(New_object)


saveRDS(New_object , file.path(PATH_ANALYSIS_OUTPUT , "Updated_Objects", "sample_NKTumor_seuratv4.rds"))
