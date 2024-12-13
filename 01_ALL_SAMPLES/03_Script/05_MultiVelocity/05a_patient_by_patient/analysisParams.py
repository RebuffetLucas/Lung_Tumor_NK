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

ANALYSIS_04a0_STEP_NAME = "04_SCENICplus_Analysis/04a0_Step0_R_Extract_Metadata"
ANALYSIS_04a1_STEP_NAME = "04_SCENICplus_Analysis/04a1_Step1_R_cisTopic_Prepro_ATACseq"
ANALYSIS_04b_STEP_NAME = "04_SCENICplus_Analysis/04b_Step2_Anndata_Prepro_scRNAseq"

ANALYSIS_STEP_NAME = "05_MultiVelocity/05a_patient_by_patient"

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
PATH_TO_THE_04d_OUPUT_FOLDER = os.path.join( PATH_EXPERIMENT_OUTPUT, "04_SCENICplus_Analysis/04d_Step4_SnakeMake_Prepare")

PATH_TO_THE_Out_01b_Neighborhood_Graph_Per_Indiv = os.path.join( PATH_ANALYSIS_OUTPUT, "Out_01b_Neighborhood_Graph_Per_Indiv")

######## CONSTANTS USED BY THE ANALYSIS STEP
######## Add below all the constants required for the analysis

CELL_TYPE_COLNAME = "Fine_Annotation_WNN"
CELL_TYPE_COLNAME2 = "Coarse_Annotation_WNN"

SAMPLE_ID_COLNAME = "sample"
ANALYSIS_SAMPLES_NAME = "NK_TUMOR"

STOP_THE_NOTEBOOK_HERE = True
DO_RUN_SCRUBLET = False

#Genes of interest after multivelocity run
GENE_LIST_OF_INTEREST = ['XCL1', 'ZEB1', 'ZEB2', 'RUNX2', 'IRF8', 'BACH2', 'ITGA1','ITGAX' , 'DNAJB1' , 'CD69','FOXO3' , 'NR4A3' ]


# FT observation
FT_OF_INTEREST= "ZEB1"
GROUP_OF_INTEREST = "NK_ENTPD1"

#Embeddings
TABLE_EMBEDDING_OF_INTEREST = "RNA_UMAP_TABLE.csv"
NUMBER_CPU = 20

#Dictionnary to rename the clusters (refer to adata.obs[CELL_TYPE_COLNAME] order)
old_to_new_mapping= {'ILC1': 'ILC1',
 'ILC3': 'ILC3',
 'NK1_CCL4': "NK1_CCL4",
 'NK1_FGFBP2': "NK1_FGFBP2",
 'NK2_XCL1': 'NK2_XCL1',
 'NK_CAMK4': "NK3_CAMK4",
 'NK_DNAJB1': "NK1_DNAJB1",
 'NK_ENTPD1': 'NK3_ENTPD1',
 'NK_GZMK': "NK2_GZMK",
 'NK_MKI67': "NK_MKI67",
 'NK_NFKB1': "NK_NFKB1"}

#Filters:
ADATA_RNA_MIN = 300
ADATA_RNA_MAX = 20000

ADATA_ATAC_MIN = 2000
ADATA_ATAC_MAX = 60000

# Define cluster names and corresponding colors
cluster_names = ["NK1_FGFBP2", "NK1_CCL4", "NK1_DNAJB1", "NK_NFKB1", "NK2_XCL1", "NK2_GZMK", "NK3_CAMK4", "NK3_ENTPD1", "NK_MKI67", "ILC1", "ILC3"]
color_codes = ["#E15759", "#FF9D9A", "#FFCC00", "#499894", "#86BCB6", "#59A14F", "#8CD17D", "#9966CC", "#CC99FF", "#6699FF", "#CCCCCC"]

#Sample of interest
SAMPLE_OF_INTEREST = "CSS1"
populations_of_interest = ["NK1_FGFBP2", "NK1_CCL4", "NK1_DNAJB1", "NK_NFKB1", "NK2_XCL1", "NK2_GZMK", "NK3_CAMK4", "NK3_ENTPD1", "NK_MKI67"] #Populations to focus on for trajectory analysis
