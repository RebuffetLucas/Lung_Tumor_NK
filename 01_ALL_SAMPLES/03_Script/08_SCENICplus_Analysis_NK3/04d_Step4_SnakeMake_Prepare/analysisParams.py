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
ANALYSIS_04a0_STEP_NAME = "08_SCENICplus_Analysis_NK3_NK3/04a0_Step0_R_Extract_Metadata"

ANALYSIS_04a1_STEP_NAME = "08_SCENICplus_Analysis_NK3/04a1_Step1_R_cisTopic_Prepro_ATACseq"

ANALYSIS_STEP_NAME = "08_SCENICplus_Analysis_NK3/04d_Step4_Anndata_Prepro_scRNAseq"


# This is the literal title of the analysis step. It will be shown at the beginning
# of the HTML report
# Example : LITERAL_TITLE = "Quality Control, normalization and clustering"
LITERAL_TITLE = "SnakeMake_Step"

# This is the path to the analysis step output folder. It will be automatically
# created at first analysis launch
PATH_ANALYSIS_OUTPUT = os.path.join( PATH_EXPERIMENT_OUTPUT, ANALYSIS_STEP_NAME)


PATH_TO_FRAGMENTS_MATRIX_CSV = os.path.join( PATH_EXPERIMENT_OUTPUT, ANALYSIS_04a1_STEP_NAME, "count_matrix.csv")
PATH_TO_CELLDATA_CSV = os.path.join( PATH_EXPERIMENT_OUTPUT, ANALYSIS_04a0_STEP_NAME, "cell_data.csv")

PATH_TO_FRAGMENTS_FILES = os.path.join( PATH_EXPERIMENT_RAWDATA, "Fragment_files")
PATH_TO_BLACK_LIST = os.path.join( PATH_EXPERIMENT_REFERENCE, "pycisTopic/blacklist/hg38-blacklist.v2.bed")

PATH_TO_INITIAL_ANNDATA_OBJECT = os.path.join( PATH_EXPERIMENT_RAWDATA, "RDS_FILES", "NK_Tumor.h5ad")
PATH_TO_THE_04a1_OUPUT_FOLDER = os.path.join( PATH_EXPERIMENT_OUTPUT, ANALYSIS_04a1_STEP_NAME)

PATH_TO_VELOCITY_FOLDER = os.path.join( PATH_EXPERIMENT_RAWDATA, "velocity", "velocity")

######## CONSTANTS USED BY THE ANALYSIS STEP
######## Add below all the constants required for the analysis


CELL_TYPE_COLNAME = "Fine_Annotation_WNN"
SAMPLE_ID_COLNAME = "sample"

ANALYSIS_SAMPLES_NAME = "NK_TUMOR"

STOP_THE_NOTEBOOK_HERE = True
DO_RUN_SCRUBLET = False


