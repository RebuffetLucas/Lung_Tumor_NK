## @knitr TransferOfLabel_StandardSeurat

#Load V2 data ref
PBMC_Meta_NK_V2 = readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/PBMC_V2_VF1_AllGenes_NewNames.rds")

if (DO_WORK_ON_11Kgenes ==TRUE){
  intersect_genes = readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/intersect_genes.rds")
  PBMC_Meta_NK_V2 = subset(PBMC_Meta_NK_V2, features= intersect_genes)
}



#print(DimPlot(PBMC_Meta_NK_V2, cols = palette))

PBMC_Meta_NK_V2$SecondClust = PBMC_Meta_NK_V2$FirstClust
levels(PBMC_Meta_NK_V2$FirstClust)= c("NK1", "NK1", "NK1", "NKint", "NK2", "NK3")
PBMC_Meta_NK_V2 = SetIdent(PBMC_Meta_NK_V2, value = "FirstClust")
#DimPlot(PBMC_Meta_NK_V2, cols = palette2)

DimPlot(PBMC_Meta_NK_V2, cols = palette2)

# Load Query
PBMC_Blood_All = readRDS(PATH_QUERY_RDS)

# Subseting Query if necessary
#Tissue subseting
if(SUBSET_META_TISSUE==TRUE){
  PBMC_Blood_subset = subset(PBMC_Blood_All, subset = meta_tissue %in%  META_TISSUE_SUBSETS)
  PBMC_Blood_subset$meta_tissue = droplevels(PBMC_Blood_subset$meta_tissue )
}else{
  PBMC_Blood_subset = PBMC_Blood_All
}

#Pick a subset of meta_histo to integrate (or not). For double criteria subsetting
if(SUBSET_META_HISTO==TRUE){
PBMC_Blood_subset = subset(PBMC_Blood_subset, idents =  META_SUBSET)
}else{
PBMC_Blood_subset = PBMC_Blood_subset
}

#Downsample randomly (or not)
if(SUBSET_DOWNSAMPLE==TRUE){
PBMC_Blood_subset = subset(PBMC_Blood_subset, downsample = N_DOWNSAMPLE)
}



#Keep genes intersection for each Dataset
genes_to_keep = intersect (rownames(PBMC_Meta_NK_V2)  , rownames(PBMC_Blood_subset ))


#Harmonize data before merging
PBMC_Blood_subset$Project = "Query"
PBMC_Blood_All@meta.data[["orig.ident"]] = PBMC_Blood_All@meta.data[[QUERY_ORIG_IDENT_COLNAME]]
PBMC_Blood_All@meta.data[["celltype"]] =  PBMC_Blood_All@meta.data[[QUERY_CLUSTER_COLNAME]]


#Merge data
data_all =  merge(PBMC_Blood_subset, y = PBMC_Meta_NK_V2, add.cell.ids = c("Query", "Meta_V2"))
data_all = subset(data_all, features= genes_to_keep)


#Normalize
data_all_list = SplitObject(data_all, split.by = DATA_SPLIT) #May be we can change to not dataset but sample for building the ref

for (i in 1:length(data_all_list)) {
  data_all_list[[i]] <- NormalizeData(data_all_list[[i]], verbose = FALSE) #Maybe we need to change to MultiBatchNorm Here
  data_all_list[[i]] <- FindVariableFeatures(data_all_list[[i]], selection.method = "vst", 
                                             nfeatures = 2000, verbose = FALSE)
}

if(DO_PREPARE_REFERENCE==TRUE){
  
#Integration
reference.list = data_all_list[REFERENCE_FOR_INTEGRATION]
Integration.anchors= FindIntegrationAnchors(object.list = reference.list, dims= 1:30)
data_integrated =  IntegrateData(anchorset = Integration.anchors, dims = 1:30)

#Dataviz of the integrated
DefaultAssay(data_integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
data_integrated <- ScaleData(data_integrated, verbose = FALSE)
data_integrated <- RunPCA(data_integrated, npcs = 30, verbose = FALSE)
data_integrated <- RunUMAP(data_integrated, reduction = "pca", dims = 1:30)

data_integrated$FirstClust = as.factor(data_integrated$FirstClust)
saveRDS(data_integrated , PATH_TO_REFERENCE )
} else {
  data_integrated = readRDS(PATH_TO_REFERENCE)
}

QUERY_DATASETS = setdiff(names(data_all_list) , REFERENCE_FOR_INTEGRATION)

DimPlot(data_integrated)


data.query = merge(data_all_list[QUERY_DATASETS][[1]] , y = data_all_list[QUERY_DATASETS][-1]) #Merge all the query together

#Projecting on the integrated reference

data.anchors <- FindTransferAnchors(reference = data_integrated, query = data.query, 
                                        dims = 1:30)
predictions <- TransferData(anchorset = data.anchors, refdata = data_integrated$FirstClust , 
                            dims = 1:30)

data.query <- AddMetaData(data.query, metadata = predictions)

saveRDS(data.query$predicted.id , file.path(PATH_ANALYSIS_OUTPUT, "predicted_id_4pops.rds") )
