#!/usr/bin/env python

#Step1a, non interactive and loop version
#Set up the env'
import os
import scipy
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import multivelo as mv
import matplotlib.pyplot as plt
import sys
import seaborn as sns
from comm import create_comm

# Determine the folder in which the code is executed
WORKING_DIR = os.getcwd()
sys.path.append(os.path.abspath( WORKING_DIR))

# Run the params codes
# Load global parameters
exec(open("../../globalParams.py").read())
# Load sample parameters
exec(open("../../sampleParams.py").read())
# Load analysis parameters
exec(open("./analysisParams.py").read())
# scvelo and panda parameters

scv.settings.verbosity = 3
scv.settings.presenter_view = True
scv.set_figure_params('scvelo')
pd.set_option('display.max_columns', 100)
pd.set_option('display.max_rows', 200)
np.set_printoptions(suppress=True)

#Looping on each of the samples:

# Define the path
# Get the list of folder names
SAMPLES_NAMES_LIST = [name for name in os.listdir(PATH_TO_VELOCITY_FOLDER) if os.path.isdir(os.path.join(PATH_TO_VELOCITY_FOLDER, name))]

for SAMPLE_OF_INTEREST in SAMPLES_NAMES_LIST :
        
    #Reading the spliced and unspliced counts
    # Read the adata object
    #adata = sc.read_h5ad("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Lung_Basel/01_ALL_SAMPLES/05_Output/04_SCENICplus_Analysis/04b_Step2_Anndata_Prepro_scRNAseq/adata.h5ad") #Need to improve this
    adata = sc.read_h5ad( os.path.join( PATH_EXPERIMENT_OUTPUT , ANALYSIS_04b_STEP_NAME , "adata.h5ad"))

    ## Update adata name
    # Create a mapping dictionary called old_to_new_mapping from the old categories to the new categories => This is directly done in analysis params.

    # Replace the old categories with the new ones
    adata.obs[CELL_TYPE_COLNAME] = adata.obs[CELL_TYPE_COLNAME].map(old_to_new_mapping)
    # Update the categories to reflect the new set
    adata.obs[CELL_TYPE_COLNAME] = adata.obs[CELL_TYPE_COLNAME].astype('category')


    ###################
    #Simplify the following code in a loop
    # load dimensional reductions:
    #RNA = pd.read_csv("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/Lung_Basel/01_ALL_SAMPLES/05_Output/04_SCENICplus_Analysis/04a0_Step0_R_Extract_Metadata/Embeddings/RNA_UMAP_TABLE.csv")
    RNA = pd.read_csv(os.path.join( PATH_EXPERIMENT_OUTPUT , ANALYSIS_04a0_STEP_NAME, "Embeddings/RNA_UMAP_TABLE.csv"), index_col = 0)
    RNA.index = [f"{cell.split('_')[1]}_{cell.split('_')[0]}" for cell in RNA.index]
    RNA = RNA.reindex(adata.obs.index)

    ATAC = pd.read_csv(os.path.join( PATH_EXPERIMENT_OUTPUT , ANALYSIS_04a0_STEP_NAME, "Embeddings/ATAC_UMAP_TABLE.csv"), index_col = 0)
    ATAC.index = [f"{cell.split('_')[1]}_{cell.split('_')[0]}" for cell in ATAC.index]
    ATAC = ATAC.reindex(adata.obs.index)

    WNN = pd.read_csv(os.path.join( PATH_EXPERIMENT_OUTPUT , ANALYSIS_04a0_STEP_NAME, "Embeddings/WNN_UMAP_TABLE.csv"), index_col = 0)
    WNN.index = [f"{cell.split('_')[1]}_{cell.split('_')[0]}" for cell in WNN.index]
    WNN = WNN.reindex(adata.obs.index)


    # set them
    adata.obsm['X_rna_umap'] = RNA.to_numpy()
    adata.obsm['X_atac_umap'] = ATAC.to_numpy()
    adata.obsm['X_wnn_umap'] = WNN.to_numpy()
    ###############################

    #Set the seed
    np.random.seed(1234)  


    
    #Load the Loom object
    # Specify the folder containing the .loom file
    folder_path = os.path.join(PATH_TO_VELOCITY_FOLDER ,SAMPLE_OF_INTEREST)

    # Find the .loom file in the folder
    loom_files = [f for f in os.listdir(folder_path) if f.endswith('.loom')]

    # Check if there is exactly one .loom file and read it
    if len(loom_files) == 1:
        loom_file_path = os.path.join(folder_path, loom_files[0])
        adata_rna = scv.read(loom_file_path, cache=True)
        print(f"Loaded {loom_files[0]}")
    else:
        print("Error: There should be exactly one .loom file in the folder.")

    # Modify the observation names
    adata_rna.obs_names = [x.split(':')[1][:-1] + '-1_' + x.split(':')[0] for x in adata_rna.obs_names]
    # Ensure variable names are unique
    adata_rna.var_names_make_unique()
    
    ### Have a quick look at where to put the threshold for total count per cell
    
    # Calculate total counts per cell
    total_counts = adata_rna.obs['n_counts'] if 'n_counts' in adata_rna.obs else adata_rna.X.sum(axis=1).A1

    # Calculate the 99th percentile to limit the x-axis
    x_limit = np.percentile(total_counts, 99)

    # Plot a histogram of total counts with more bins
    plt.figure(figsize=(10, 6))
    sns.histplot(total_counts, bins=100, kde=True)  # Increased bins for more granularity
    plt.xlabel('Total Counts per Cell')
    plt.ylabel('Number of Cells')
    plt.title('Distribution of Total Counts per Cell')
    plt.axvline(x=ADATA_RNA_MIN, color='red', linestyle='--', label=f"min_counts={ADATA_RNA_MIN}")
    plt.axvline(x=ADATA_RNA_MAX, color='green', linestyle='--', label= f"max_counts={ADATA_RNA_MAX}")
    plt.xlim(0, x_limit)  # Limit the x-axis to the 99th percentile
    plt.legend()
    plt.show()
    
    #Apply filter on cells
    sc.pp.filter_cells(adata_rna, min_counts=ADATA_RNA_MIN)
    sc.pp.filter_cells(adata_rna, max_counts=ADATA_RNA_MAX)
    
    #Keep only the top 1000 genes
    # Top 1000 variable genes are used for downstream analyses.
    scv.pp.filter_and_normalize(adata_rna, min_shared_counts=10, n_top_genes=1000)
    
    
    #Compare filtered loom object with the total object
    # Find the intersection of observation names
    common_obs_names = set(adata.obs_names).intersection(adata_rna.obs_names)

    # Convert to a list if needed
    common_obs_names = list(common_obs_names)

    # Output the number of common observation names
    print(f"Number of common observation names: {len(common_obs_names)}")

    # Optionally, display the common observation names
    #print(common_obs_names)

    # Filter adata.obs for rows where the sample is "sample of interest and get the row names
    sample_cells = adata.obs_names[adata.obs['sample'] == SAMPLE_OF_INTEREST]

    # Convert to a list if needed
    sample_cells = list(sample_cells)

    # Output the row names
    print(f"Number of cells in the same patient in adata object: {len(sample_cells)}")
    
    
    ### Keep only the data that are present in the final object and insert the metadata in adata_rna
    # Step 1: Subset adata_rna to only include common observation names
    adata_rna = adata_rna[adata_rna.obs_names.isin(common_obs_names)].copy()

    # Step 2: Create adata.obs2 containing only common observation names from adata.obs
    adata_obs2 = adata.obs.loc[adata.obs_names.isin(common_obs_names)].copy()

    # Step 3: Reorder adata_obs2 to match the order in adata_rna
    adata_obs2 = adata_obs2.reindex(adata_rna.obs_names)

    # Step 4: Add the information from adata_obs2 to adata_rna.obs
    adata_rna.obs = pd.concat([adata_rna.obs, adata_obs2], axis=1)
    
    #Focusing on pops of interest to elaborate the trajectories
    adata_rna = adata_rna[adata_rna.obs[CELL_TYPE_COLNAME].isin(populations_of_interest)]
    
    
    
    #Preprocessing the ATAC part
    #Load the atac object
    # Specify the folder containing the .loom file
    folder_path = os.path.join(
        PATH_TO_VELOCITY_FOLDER,
        SAMPLE_OF_INTEREST,
        f"{SAMPLE_OF_INTEREST}_filtered_feature_bc_matrix"
    )

    adata_atac = sc.read_10x_mtx(folder_path, var_names='gene_symbols', cache=True, gex_only=False)
    adata_atac = adata_atac[:,adata_atac.var['feature_types'] == "Peaks"]
    
    # We aggregate peaks around each gene as well as those that have high correlations with promoter peak or gene expression.
    # Peak annotation contains the metadata for all peaks.
    # Feature linkage contains pairs of correlated genomic features.
    adata_atac = mv.aggregate_peaks_10x(adata_atac,
                                        os.path.join(   PATH_TO_VELOCITY_FOLDER,    SAMPLE_OF_INTEREST,   f"{SAMPLE_OF_INTEREST}_atac_peak_annotation.tsv"),
                                        os.path.join(   PATH_TO_VELOCITY_FOLDER,    SAMPLE_OF_INTEREST,   f"{SAMPLE_OF_INTEREST}_feature_linkage.bedpe")
    )
                        

    # Let's examine the total count distribution and remove outliers.
    plt.hist(adata_atac.X.sum(1), bins=100, range=(0, ADATA_ATAC_MAX));
    
    
    sc.pp.filter_cells(adata_atac, min_counts=ADATA_ATAC_MIN)
    sc.pp.filter_cells(adata_atac, max_counts=ADATA_ATAC_MAX)
    
    # We normalize aggregated peaks with TF-IDF.
    mv.tfidf_norm(adata_atac)
    
    ## Finding shared barcodes and features between RNA and ATAC
    # Update obs_names by adding "_SAMPLE_OF_INTEREST" at the end of each barcode
    adata_atac.obs_names = [f"{barcode}_{SAMPLE_OF_INTEREST}" for barcode in adata_atac.obs_names]

    # Verify the changes
    #print(adata_atac.obs_names[:5])  # Display the first few updated names for confirmation

    shared_cells = pd.Index(np.intersect1d(adata_rna.obs_names, adata_atac.obs_names))
    shared_genes = pd.Index(np.intersect1d(adata_rna.var_names, adata_atac.var_names))
    len(shared_cells), len(shared_genes)
    
    
    #Reload the data in and carry on with a subset of cells
    # Reloading the adata_rna and adding its info
    #Load the Loom object
    # Specify the folder containing the .loom file
    folder_path = os.path.join(PATH_TO_VELOCITY_FOLDER ,SAMPLE_OF_INTEREST)

    # Find the .loom file in the folder
    loom_files = [f for f in os.listdir(folder_path) if f.endswith('.loom')]

    # Check if there is exactly one .loom file and read it
    if len(loom_files) == 1:
        loom_file_path = os.path.join(folder_path, loom_files[0])
        adata_rna = scv.read(loom_file_path, cache=True)
        print(f"Loaded {loom_files[0]}")
    else:
        print("Error: There should be exactly one .loom file in the folder.")

    # Modify the observation names
    adata_rna.obs_names = [x.split(':')[1][:-1] + '-1_' + x.split(':')[0] for x in adata_rna.obs_names]
    # Ensure variable names are unique
    adata_rna.var_names_make_unique()

    ## Adding the metadata in
    # Find the intersection of observation names
    common_obs_names = set(adata.obs_names).intersection(adata_rna.obs_names)

    # Convert to a list if needed
    common_obs_names = list(common_obs_names)

    # Output the number of common observation names
    print(f"Number of common observation names: {len(common_obs_names)}")

    # Optionally, display the common observation names
    #print(common_obs_names)

    # Filter adata.obs for rows where the sample is "sample of interest and get the row names
    sample_cells = adata.obs_names[adata.obs['sample'] == SAMPLE_OF_INTEREST]

    # Convert to a list if needed
    sample_cells = list(sample_cells)

    # Output the row names
    print(f"Number of cells in the same patient in adata object: {len(sample_cells)}")

    # Step 1: Subset adata_rna to only include common observation names
    adata_rna = adata_rna[adata_rna.obs_names.isin(common_obs_names)].copy()

    # Step 2: Create adata.obs2 containing only common observation names from adata.obs
    adata_obs2 = adata.obs.loc[adata.obs_names.isin(common_obs_names)].copy()

    # Step 3: Reorder adata_obs2 to match the order in adata_rna
    adata_obs2 = adata_obs2.reindex(adata_rna.obs_names)

    # Step 4: Add the information from adata_obs2 to adata_rna.obs
    adata_rna.obs = pd.concat([adata_rna.obs, adata_obs2], axis=1)
    
    #Subsetting the two objects to keep only the cells and genes of interest
    adata_rna = adata_rna[shared_cells, shared_genes]
    adata_atac = adata_atac[shared_cells, shared_genes]
    
    
    #Running standard procedure
    scv.pp.normalize_per_cell(adata_rna)
    scv.pp.log1p(adata_rna)
    scv.pp.moments(adata_rna, n_pcs=30, n_neighbors=50)
    
    #Make sure that data annotation is categorial
    adata_rna.obs[CELL_TYPE_COLNAME] = adata_rna.obs[CELL_TYPE_COLNAME].astype('category')
    
    scv.tl.umap(adata_rna)
    scv.pl.umap(adata_rna, color= CELL_TYPE_COLNAME)
    
    import scvelo as scv

    # Ensure that the color map is defined for the specific clusters
    cluster_color_map = dict(zip(cluster_names, color_codes))

    # Optionally, add titles and labels for better visualization
    scv.pl.umap(
        adata_rna, 
        color=CELL_TYPE_COLNAME, 
        palette=cluster_color_map, 
        size=50, 
        legend_loc='right', 
        legend_fontsize=10,
        title='UMAP of Cell Types', 
        xlabel='UMAP 1', 
        ylabel='UMAP 2', 
        frameon=False
    )

    #Inject WNN UMAP dimensions in the object
    WNN = pd.read_csv(os.path.join( PATH_EXPERIMENT_OUTPUT , ANALYSIS_04a0_STEP_NAME, "Embeddings/WNN_UMAP_TABLE.csv"), index_col = 0)
    # Modify the row names to the desired format
    WNN.index = [f"{cell.split('_')[1]}_{cell.split('_')[0]}" for cell in WNN.index]
    WNN = WNN.reindex(adata_rna.obs.index)

    #Inject into the object
    # set them in the umap slot
    adata_rna.obsm['X_umap'] = WNN.to_numpy()
    
    # Plot the UMAP using the 'X_wnn_umap' embedding
    # Optionally, add titles and labels for better visualization
    scv.pl.umap(
        adata_rna, 
        color=CELL_TYPE_COLNAME, 
        palette=cluster_color_map, 
        size=50, 
        legend_loc='right', 
        legend_fontsize=10,
        title='WNN UMAP of Cell Types', 
        xlabel='UMAP 1', 
        ylabel='UMAP 2', 
        frameon=False
    )

    #Saving the filtered cells
    # Write out filtered cells and prepare to run Seurat WNN --> R script can be found on Github.

    # Construct the file path
    file_path = os.path.join(
        PATH_ANALYSIS_OUTPUT, 
        "Out_01a_Filteredcells_per_indiv_sample", "filtered_cells",
        f"{SAMPLE_OF_INTEREST}_filtered_cells.txt"
    )
    
    # Create directories if they don't exist
    os.makedirs(os.path.dirname(file_path), exist_ok=True)

    # Save the obs_names to the specified file path
    adata_rna.obs_names.to_frame().to_csv(file_path, header=False, index=False)

    # Print the file path to confirm where the file was saved
    print(f"Filtered cells saved to: {file_path}")
    
    #save the adata_rna and adata_atac objects
    # Construct the file path
    file_path = os.path.join(
        PATH_ANALYSIS_OUTPUT, 
        "Out_01a_Filteredcells_per_indiv_sample", "adata_rna",
        f"{SAMPLE_OF_INTEREST}_adata_rna.h5ad"
    )
    
    # Create directories if they don't exist
    os.makedirs(os.path.dirname(file_path), exist_ok=True)

    # Save the obs_names to the specified file path
    adata_rna.write(file_path)

    #save the adata_rna and adata_atac objects
    # Construct the file path
    file_path = os.path.join(
        PATH_ANALYSIS_OUTPUT, 
        "Out_01a_Filteredcells_per_indiv_sample", "adata_atac",
        f"{SAMPLE_OF_INTEREST}_adata_atac.h5ad"
    )
    
    # Create directories if they don't exist
    os.makedirs(os.path.dirname(file_path), exist_ok=True)

    # Save the obs_names to the specified file path
    adata_atac.write(file_path)