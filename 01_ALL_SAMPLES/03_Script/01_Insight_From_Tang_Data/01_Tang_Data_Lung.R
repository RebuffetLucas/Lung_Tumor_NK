

PBMC_All = LoadH5Seurat("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/01_Reference/Tang_Data/comb_CD56_CD16_NK.h5seurat")

PBMC_Lung = subset( PBMC_All , subset = meta_histology== "Lung Cancer(LC)")
PBMC_Lung@meta.data = droplevels(PBMC_Lung@meta.data)
