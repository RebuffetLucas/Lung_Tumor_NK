###############################################################################
# This file defines ANALYSIS parameters as global variables that will be loaded
# before analysis starts. It should define common parameters used by the current
# analysis


ANALYSIS_STEP_NAME = "06_Cytosig"
PATH_ANALYSIS_OUTPUT = file.path( PATH_EXPERIMENT_OUTPUT, ANALYSIS_STEP_NAME)

LITERAL_TITLE = "NK_Lung_Tumor"


MIN.CELLS=3
MIN.FEATURES=200
RIBO_THRESHOLD=10
MITO_THRESHOLD=0
N_GENES_VARIABLES=2000
DO_SCALE=FALSE
DO_CENTER=TRUE
DIM_PCA=50
DIM_UMAP=50
RESOLUTION=0.6
FDRCUTOFF=0.05
FIND_ALL_MARKERS_LOGFC=0.25


DOWNSAMPLE_ATAC_FOR_TEST = TRUE
DOWN_SAMPLE_PEAKS = 5000
DOWN_SAMPLE_CELLS = 5000


# Variables specific for cytosig analysis
#celltype_col = 
#value_col = 

numOfCellsPerRun <- 10000

cytokines <- c('IL15','IL2','IL21','IL12','IL18','IL10','IL27','IFN1','IFNG','TGFB1','TNFA','PGE2','TRAIL')





