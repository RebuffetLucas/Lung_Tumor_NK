###############################################################################
# This file defines ANALYSIS parameters as global variables that will be loaded
# before analysis starts. It should define common parameters used by the current
# analysis



ANALYSIS_STEP_NAME = "12_Trajectory_Characterization"
PATH_ANALYSIS_OUTPUT = file.path( PATH_EXPERIMENT_OUTPUT, ANALYSIS_STEP_NAME)

LITERAL_TITLE = "NK_Lung_Cancer"

SAMPLE_NAME = "SAMPLE_NK3_Only_Preprocessed"


# Monocle params
CLUSTERS_TO_REMOVE = c("NK3","NK2")

LIST_PARAM_CONTROL = list(ncenter=165)



root_pr_nodes = "Y_14"

#For graph test:
#Q_Value_Limit = 0.00005
#Coefficients
#k=9

NUMBER_GENES_DYNAMIC_HEATMAP = 100
NUMBER_GENES_DYNAMIC_HEATMAP_SMALL = 50

# Others
MIN.CELLS=3
MIN.FEATURES=200
RIBO_THRESHOLD=10
MITO_THRESHOLD=0
N_GENES_VARIABLES=2000
DO_SCALE=FALSE
DO_CENTER=TRUE
DIM_PCA=40
DIM_LSI=40
DIM_UMAP=40
RESOLUTION=0.6
FDRCUTOFF=0.05
FIND_ALL_MARKERS_LOGFC=0.25


DOWNSAMPLE_ATAC_FOR_TEST = FALSE
DOWN_SAMPLE_PEAKS = 20000
DOWN_SAMPLE_CELLS = 20000

SEED= 123

NCORES_LDA_MODELS = 10


DO_EXTRACT_EMBED = FALSE

DO_EXTRACT_NEIGBORS = FALSE

# For Preprocessing step
variable_to_sort = "Fine_Annotation_WNN"
Clusters_To_Extract = c("NK_GZMK", "NK_CAMK4",   "NK_ENTPD1")

# Dictionary mapping samples to orig.ident
sample_to_orig <- list(
  CSS1 = "BS-840",
  CSS10 = "BS-897",
  CSS13 = "BS-1175",
  CSS16 = "BS-1198",
  CSS19 = "BS-1308",
  CSS21 = "BS-1314",
  CSS23 = "BS-1319",
  CSS25 = "BS-1318",
  CSS27 = "BS-1322",
  CSS4 = "BS-824",
  CSS7 = "BS-889"
)


# Define cluster names and corresponding colors
cluster_names = c("NK3_GZMK", "NK3_CAMK4", "NK3_ENTPD1")
color_codes = c( "#59A14F", "#8CD17D", "#9966CC")

col_firstCl=c(NK3_GZMK = "#59A14F", 
              NK3_CAMK4 = "#8CD17D",
              NK3_ENTPD1 = "#9966CC")

Annotation_of_interest = "Fine_Annotation_WNN_2"


Choose_root_method = "Supervised" # "Automatic"

Choose_feature_gene_activty_method = "q95" # "Top_VST" or "q95"



#For regulons:
REGULON_TO_USE = "Filtered_Activators" #"Filtered_Repressors" # "All" , "Filtered_Activators"

Keep_only_plus_plus_regulons = TRUE







