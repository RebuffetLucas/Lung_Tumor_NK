#This code aims at extracting a series of reference signatures that will be used in the project
#Sce object with all the cells
PBMCscaled= readRDS(file= "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK3/00_RawData/Seurat_Obj_other_samples/Sent_By_Janine/merged.rds")

#Convert it to seurat Object
colnames(PBMCscaled)=make.unique(colnames(PBMCscaled))
PBMCscaled= as.Seurat(PBMCscaled)

#Composition before downsampling
table(PBMCscaled@meta.data[["orig.ident"]], PBMCscaled@meta.data[["Chemistry"]])

#Subseting the data
PBMCscaled2  =subset(PBMCscaled, subset= Chemistry == "V2" & Dataset == "Dataset4"  & nkg2c == "NKG2Cpos")

table(PBMCscaled2@meta.data[["orig.ident"]], PBMCscaled2@meta.data[["Chemistry"]])


#Setident to orig.ident:
PBMCscaled2 = SetIdent(PBMCscaled2, value = "orig.ident")

#Have a look at the cells
DimPlot(PBMCscaled2)

PBMCscaled2= NormalizeData(PBMCscaled2)
#Find Markers CMVpos vs CMVneg:
CMVmarkers = FindMarkers(PBMCscaled2, ident.1 = c("CMVpos2_donorC" , "CMVpos3_donorD" ,  "CMVpos4_donorE") , ident.2 = c("CMVneg1_donorA", "CMVneg2_donorB") )

saveRDS(CMVmarkers, file.path(PATH_EXPERIMENT_OUTPUT, "02_Preparing_DataTables_For_Scoring","Signature_HCMVpos_NKG2Cpos_vsHCMVneg_NKG2Cpos.rds"))
