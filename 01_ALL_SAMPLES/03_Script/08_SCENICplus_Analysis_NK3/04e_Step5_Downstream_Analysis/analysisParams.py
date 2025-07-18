###############################################################################
# This file defines ANALYSIS parameters as global variables that will be loaded
# before analysis starts. It should define common parameters used by the current
# analysis

import os

######## CONSTANTS ON THE ANALYSIS STEP IDENTITY
######## Modify the ANALYSIS_STEP_NAME and LITERAL_TITLE variable

# This is the name of the analysis step. Use the same name as the folder
# name of the analysis step
# Example : ANALYSIS_STEP_NAME = "02_GlobalHeterogeneity"

ANALYSIS_04a0_STEP_NAME = "08_SCENICplus_Analysis_NK3/04a0_Step0_R_Extract_Metadata"
ANALYSIS_04a1_STEP_NAME = "08_SCENICplus_Analysis_NK3/04a1_Step1_R_cisTopic_Prepro_ATACseq"



ANALYSIS_STEP_NAME = "08_SCENICplus_Analysis_NK3/04e_Step5_Downstream_Analysis"

# This is the literal title of the analysis step. It will be shown at the beginning
# of the HTML report
# Example : LITERAL_TITLE = "Quality Control, normalization and clustering"
LITERAL_TITLE = "Downstream_Analysis_SCenicplus_Outputs"

# This is the path to the analysis step output folder. It will be automatically
# created at first analysis launch
PATH_ANALYSIS_OUTPUT = os.path.join( PATH_EXPERIMENT_OUTPUT, ANALYSIS_STEP_NAME)


PATH_TO_FRAGMENTS_MATRIX_CSV = os.path.join( PATH_EXPERIMENT_OUTPUT, ANALYSIS_04a1_STEP_NAME, "count_matrix.csv")
PATH_TO_CELLDATA_CSV = os.path.join( PATH_EXPERIMENT_OUTPUT, ANALYSIS_04a0_STEP_NAME, "cell_data.csv")

PATH_TO_FRAGMENTS_FILES = os.path.join( PATH_EXPERIMENT_RAWDATA, "Fragment_files")
PATH_TO_BLACK_LIST = os.path.join( PATH_EXPERIMENT_REFERENCE, "pycisTopic/blacklist/hg38-blacklist.v2.bed")

PATH_TO_INITIAL_ANNDATA_OBJECT = os.path.join( PATH_EXPERIMENT_RAWDATA, "RDS_FILES", "NK_Tumor.h5ad")

PATH_TO_VELOCITY_FOLDER = os.path.join( PATH_EXPERIMENT_RAWDATA, "velocity", "velocity")

PATH_TO_THE_04d_OUPUT_FOLDER = os.path.join( PATH_EXPERIMENT_OUTPUT, "08_SCENICplus_Analysis_NK3/04d_Step4_SnakeMake_Prepare")



######## CONSTANTS USED BY THE ANALYSIS STEP
######## Add below all the constants required for the analysis


CELL_TYPE_COLNAME = "Fine_Annotation_WNN"

CELL_TYPE_COLNAME2 = "Coarse_Annotation_WNN"


SAMPLE_ID_COLNAME = "sample"

ANALYSIS_SAMPLES_NAME = "NK_TUMOR"

STOP_THE_NOTEBOOK_HERE = True
DO_RUN_SCRUBLET = False

FT_OF_INTEREST= "IRF4"

GROUP_OF_INTEREST = "NK_ENTPD1"

TABLE_EMBEDDING_OF_INTEREST = "WNN_UMAP_TABLE.csv"

NUMBER_CPU = 20

# Define cluster names and corresponding colors
cluster_names = ['NK1_FGFBP2', 'NK1_CCL4','NK_DNAJB1', 'NK2_NFKB1',   'NK2_XCL1','NK3_GZMK','NK3_CAMK4', 'NK_ENTPD1',  'NK_Prolif','ILC1', 'ILC3'] 
color_codes = ["#E15759", "#FF9D9A", "#FFCC00", "#499894", "#86BCB6", "#59A14F", "#8CD17D", "#9966CC", "#CC99FF", "#6699FF", "#CCCCCC"]

GROUP_VARIABLE_ORDER = ['ILC1', 'ILC3', 'NK_Prolif','NK1_FGFBP2', 'NK1_CCL4', 'NK2_NFKB1', 'NK2_XCL1','NK3_GZMK','NK3_CAMK4', 'NK_ENTPD1'  ]

#Filtering or not for the regulons +/+ and +/-
DO_FILTER_REGULON = True