###############################################################################
# This file defines ANALYSIS parameters as global variables that will be loaded
# before analysis starts. It should define common parameters used by the current
# analysis

######## CONSTANTS ON THE ANALYSIS STEP IDENTITY
######## Modify the ANALYSIS_STEP_NAME and LITERAL_TITLE variable

# This is the name of the analysis step. Use the same name as the folder
# name of the analysis step
# Example : ANALYSIS_STEP_NAME = "02_GlobalHeterogeneity"
ANALYSIS_STEP_NAME = "09_Label_Transfer"
SUB_ANALYSIS_STEP_NAME= "MetaNK_Label_Transfer_4pops"

# This is the literal title of the analysis step. It will be shown at the beginning
# of the HTML report
# Example : LITERAL_TITLE = "Quality Control, normalization and clustering"
LITERAL_TITLE = "Projecting Meta-NK classification in NK tumor samples"

# This is the path to the analysis step output folder. It will be automatically
# created at first analysis launch
PATH_ANALYSIS_OUTPUT = file.path( PATH_EXPERIMENT_OUTPUT, ANALYSIS_STEP_NAME, SUB_ANALYSIS_STEP_NAME)
PATH_ANALYSIS_OUTPUT_FIGURES = file.path( PATH_ANALYSIS_OUTPUT , "Figures") 
#Create the path output if it does not already exist
# Specify the directory path
dir_path <-PATH_ANALYSIS_OUTPUT
# Check if the directory exists; if not, create it
if (!file.exists(dir_path)) {
  dir.create(dir_path, recursive = TRUE)
  cat("Directory created at:", dir_path, "\n")
} else {
  cat("Directory already exists:", dir_path, "\n")
}

dir_path <-PATH_ANALYSIS_OUTPUT_FIGURES
# Check if the directory exists; if not, create it
if (!file.exists(dir_path)) {
  dir.create(dir_path, recursive = TRUE)
  cat("Directory created at:", dir_path, "\n")
} else {
  cat("Directory already exists:", dir_path, "\n")
}


######## CONSTANTS USED BY THE ANALYSIS STEP
######## Add here all the constants required for the analysis
DO_WORK_ON_11Kgenes = TRUE



#Query Dataset
PATH_QUERY_RDS = file.path(PATH_ANALYSIS_OUTPUT , "Updated_Objects", "sample_NKTumor_seuratv4.rds" ) #Needs to be updatedd
#Analysis params
DO_SUBSET= FALSE
VARIABLE_TO_SUBSET = "orig.ident"

SUBSET_META_TISSUE= FALSE #Whether to subset on some specific meta_histology categories
META_TISSUE_SUBSETS = c("All")

SUBSET_META_HISTO = TRUE #Whether to subset on some specific meta_histology categories
META_SUBSET =  c("NK1_FGFBP2" ,"NK1_CCL4"  , "NK_NFKB1" ,  "NK2"  , "NK3_GZMK",    "NK3_CAMK4" ,  "NK3_ENTPD1" , "NK_MKI67",  "NK_DNAJB1" )

SUBSET_DOWNSAMPLE = FALSE #Whether downsample Tang Data
N_DOWNSAMPLE = 2000 #Number of cells to downasample randomly


PATH_TO_REFERENCE = PATH_TO_REFERENCE = file.path("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Lyon_NK_Cytokine_Project/01_SAMPLES_ONE_BY_ONE/05_Output/02_Label_Transfer/MetaNK_Label_Transfer_4pops/reference_4pops.rds")
DO_PREPARE_REFERENCE = TRUE

#Harmonyzing the metadata before merging
QUERY_ORIG_IDENT_COLNAME = "orig.ident"
QUERY_CLUSTER_COLNAME = "Fine_Annotation_WNN"

#Building the reference
DATA_SPLIT =  "orig.ident" #Data to split by for normalization and label transfer procedure #BatchCor = samples for MetaNK and dataset for Tang



DO_REMOVE_VERY_LITTLE_IDENT_SUBSET = FALSE #If there are less than N_TRESHOLD_CELL for Diagnostic part
TRESHOLD_NEW_ID_REMOVING = 10 # Minimal number of cells identified to keep the cluster


#REFERENCE_FOR_INTEGRATION =  c("Dataset1", "Dataset2", "Dataset3", "Dataset4")

REFERENCE_FOR_INTEGRATION =  c("CMVneg1_donorA", "CMVneg2_donorB", "CMVpos2_donorC", "CMVpos3_donorD", "CMVpos4_donorE", "GSM3377678" ,"GSM3738542", "GSM3738543", "GSM5584154"  ,"GSM5584155", "GSM5584156_1", "GSM5584156_2", "GSM5584156_3")  

                               
#Defining variable features
#VARIABLE_FEATURES_METHOD = "FORCE_VF1" #FORCE_VF1 = take the intersection of VF1 and rownames as variable features
#"VARIABLE_FEATURE_BASIC" = Take classical 2000 VF within select integrations
VARIABLE_FEATURES_METHOD  = "VARIABLE_FEATURE_INTEGRATION_FEATURES" #=  Select Integration features 
VARIABLE_FEATURE_VAR_TO_SPLIT_TO_FIND_FEATURES= "orig.ident"


HARMO_USE_REF = TRUE #True or false , WHether to use a reference during Harmony Integration
HARMO_REF_TO_USE  = c("CMVneg1_donorA", "CMVneg2_donorB" ,"CMVpos2_donorC", "CMVpos3_donorD", "CMVpos4_donorE" ,"GSM3377678" ,    "GSM3738542",
                      "GSM3738543" ,    "GSM5584154" ,    "GSM5584155"  ,   "GSM5584156_1" ,  "GSM5584156_2",   "GSM5584156_3"  )


#Save RDS file?
SAVE_RDS = TRUE
RDS_DIR_OUTPUT = file.path(PATH_ANALYSIS_OUTPUT, "Updated_Object_With_predicted_id_4pops.rds")


NBCORES = 8


#Module scores
NUMBER_TOP_SCORING = 20


#Colors
palette<-c('NK1C'='#F8766D','NKint'='#8494FF',
           'NK1A'='#0CB702',
           'NK1B'='#00BFC4','NK2'='#ED68ED',
           'NK3'='#ABA300')

palette2<-c('NK1'='#F8766D','NKint'='#8494FF',
           'NK2'='#ED68ED',
           'NK3'='#ABA300')



# Seed for pseudo-random numbers
SEED = 42;


# Scaling parameters (see Seurat::ScaleData())
DATA_CENTER       = TRUE
DATA_SCALE        = TRUE
DATA_VARS_REGRESS = NULL  # c("nCount_RNA") for UMIs (NULL to ignore)

#### Analysis parameters

# Maximum number of variable features to keep
VARIABLE_FEATURES_MAXNB   = 2000;  # For analysis (PCA)
VARIABLE_FEATURES_SHOWTOP = 200;   # For table in report


# Nearest-neighbor graph construction
FINDNEIGHBORS_K = 30

# Cluster identification parameters
FINDCLUSTERS_RESOLUTION     = 0.5;
FINDCLUSTERS_USE_PCA_NBDIMS = 30;  # Number of dimensions to use from PCA results
FINDCLUSTERS_ALGORITHM      = 1;   # 1 = Louvain; 2 = Louvain with multilevel refinement; 3 = SLM; 4 = Leiden
FINDCLUSTERS_DIMS = 1:30

# Parameters for identification of marker annotations for clusters (see Seurat::FindAllMarkers())
FINDMARKERS_METHOD    = "wilcox"  # Method used to identify markers
FINDMARKERS_ONLYPOS   = TRUE     # Only consider overexpressed annotations for markers ? (if FALSE downregulated genes can also be markers)

FINDMARKERS_MINPCT    = 0.2      # Only test genes that are detected in a minimum fraction of cells in either of the two populations. Speed up the function by not testing genes that are very infrequently expressed. Default is '0.1'.
FINDMARKERS_LOGFC_THR = 0.25     # Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells. Default is '0.25'. Increasing logfc.threshold speeds up the function, but can miss weaker signals.
FINDMARKERS_PVAL_THR  = 0.001    # PValue threshold for identification of significative markers
FINDMARKERS_SHOWTOP   = 10      # Number of marker genes to show in report and tables (NULL for all)

FDRCUTOFF = 0.05

