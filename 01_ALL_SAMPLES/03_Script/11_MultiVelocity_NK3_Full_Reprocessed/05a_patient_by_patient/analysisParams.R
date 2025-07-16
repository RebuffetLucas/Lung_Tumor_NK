###############################################################################
# This file defines ANALYSIS parameters as global variables that will be loaded
# before analysis starts. It should define common parameters used by the current
# analysis


ANALYSIS_STEP_NAME = "11_MultiVelocity_NK3_Full_Reprocessed/05a_patient_by_patient"
PATH_ANALYSIS_OUTPUT = file.path( PATH_EXPERIMENT_OUTPUT, ANALYSIS_STEP_NAME)

LITERAL_TITLE = "NK_Lung_Cancer"

SAMPLE_NAME = "Intra_Tum_NK3_Only"

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


NCORES_LDA_MODELS = 10


DO_EXTRACT_EMBED = TRUE

DO_EXTRACT_NEIGBORS = TRUE

# For Preprocessing step
variable_to_sort = "Fine_Annotation_WNN_2"
Clusters_To_Extract = c("NK_GZMK", "NK_CAMK4",   "NK_ENTPD1") #To update here

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



