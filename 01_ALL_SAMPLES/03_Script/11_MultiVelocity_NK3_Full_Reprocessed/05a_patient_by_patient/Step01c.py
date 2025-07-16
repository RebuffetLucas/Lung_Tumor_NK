#!/usr/bin/env python

#Step1c, non interactive and loop version
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

#Remove samples that do not contain enought cells
if DO_EXCEPT_SOME_PATIENTS:
    SAMPLES_NAMES_LIST = [sample for sample in SAMPLES_NAMES_LIST if sample not in SAMPLE_NOT_TO_DO]

    
for SAMPLE_OF_INTEREST in SAMPLES_NAMES_LIST :
    print(SAMPLE_OF_INTEREST)
    #Read in adata_rna and adata_atac
    #adata_rna
    file_path = os.path.join(
    PATH_ANALYSIS_OUTPUT, 
    "Out_01a_Filteredcells_per_indiv_sample", "adata_rna",
    f"{SAMPLE_OF_INTEREST}_adata_rna.h5ad"
    )
    adata_rna = sc.read_h5ad(file_path)
    
    # adata_atac
    # Construct the file path
    file_path = os.path.join(
        PATH_ANALYSIS_OUTPUT, 
        "Out_01a_Filteredcells_per_indiv_sample", "adata_atac",
        f"{SAMPLE_OF_INTEREST}_adata_atac.h5ad"
    )

    # Save the obs_names to the specified file path
    adata_atac = sc.read_h5ad(file_path)
    
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
    
    
    ## Smoothing gene aggregated peaks by neighbors
    
    # Read in Seurat WNN neighbors.
    nn_idx = np.loadtxt(os.path.join( PATH_TO_THE_Out_01b_Neighborhood_Graph_Per_Indiv, SAMPLE_OF_INTEREST,"nn_idx.txt"), delimiter=',')
    nn_dist = np.loadtxt(os.path.join(PATH_TO_THE_Out_01b_Neighborhood_Graph_Per_Indiv, SAMPLE_OF_INTEREST,"nn_dist.txt"), delimiter=',')
    nn_cells = pd.Index(pd.read_csv(os.path.join( PATH_TO_THE_Out_01b_Neighborhood_Graph_Per_Indiv, SAMPLE_OF_INTEREST,"nn_cells.txt"), header=None)[0])

    # Make sure cell names match.
    np.all(nn_cells == adata_atac.obs_names)
    
    mv.knn_smooth_chrom(adata_atac, nn_idx, nn_dist)
    
    # Running multi-omic dynamical model
    #This will take a while. Parallelization is high recommended.
    mv.settings.VERBOSITY = 0

    adata_result = mv.recover_dynamics_chrom(adata_rna,
                                             adata_atac,
                                             max_iter=5,
                                             init_mode="invert",
                                             parallel=True,
                                             save_plot=False,
                                             rna_only=False,
                                             fit=True,
                                             n_anchors=500,
                                             extra_color_key=CELL_TYPE_COLNAME)
    
    
    # Save the result for use later on
    # Construct the file path
    file_path = os.path.join(
        PATH_ANALYSIS_OUTPUT,"Out_01c_MultiVeloAnalysis" ,"multi_velo_results",
        f"{SAMPLE_OF_INTEREST}_multivelo_result.h5ad"
    )

    # Save the AnnData object to the specified file path
    adata_result.write(file_path)

    # Print the file path to confirm the save location
    print(f"Result saved to: {file_path}")
    
    # Construct the file path for saving the plot
    file_path = os.path.join(
        PATH_ANALYSIS_OUTPUT, 
        "Figures/Quality_metrics_Plots", 
        f"{SAMPLE_OF_INTEREST}_Plot1.pdf"
    )

    # Create the directories if they don't exist
    os.makedirs(os.path.dirname(file_path), exist_ok=True)

    mv.pie_summary(adata_result)

    # Generate the pie summary plot
    mv.pie_summary(adata_result)

    # Save the current figure to the specified file path
    plt.savefig(file_path, format='pdf')

    # Close the figure to free memory
    plt.close()
    
    # Construct the file path for saving the plot
    file_path = os.path.join(
        PATH_ANALYSIS_OUTPUT, 
        "Figures/Quality_metrics_Plots", 
        f"{SAMPLE_OF_INTEREST}_Plot2.pdf"
    )

    # Create the directories if they don't exist
    os.makedirs(os.path.dirname(file_path), exist_ok=True)


    # Generate the switch time summary plot
    mv.switch_time_summary(adata_result)
    plt.xticks(rotation=45)  # Rotate the labels by 45 degrees (adjust the angle as needed)

    # Rotate the x-axis labels
    mv.switch_time_summary(adata_result)
    plt.xticks(rotation=45)  # Rotate the labels by 45 degrees (adjust the angle as needed)
    # Save the current figure to the specified file path
    plt.savefig(file_path, format='pdf')

    # Close the figure to free memory
    plt.close()
    
    # Construct the file path for saving the plot
    file_path = os.path.join(
        PATH_ANALYSIS_OUTPUT, 
        "Figures/Quality_metrics_Plots", 
        f"{SAMPLE_OF_INTEREST}_Plot3.pdf"
    )

    # Create the directories if they don't exist
    os.makedirs(os.path.dirname(file_path), exist_ok=True)

    mv.likelihood_plot(adata_result)

    # Generate the pie summary plot
    mv.likelihood_plot(adata_result)

    # Save the current figure to the specified file path
    plt.savefig(file_path, format='pdf')

    # Close the figure to free memory
    plt.close()
    
    
    #Computing Velocity stream and latent time
    mv.velocity_graph(adata_result)
    mv.latent_time(adata_result)
    
    # Construct the file path for saving the plot
    file_path = os.path.join(
        PATH_ANALYSIS_OUTPUT, 
        "Figures/Velocity_streams_per_individuals", 
        f"{SAMPLE_OF_INTEREST}_Velostream.svg"
    )

    # Create the directories if they don't exist
    os.makedirs(os.path.dirname(file_path), exist_ok=True)
    
    # Generate the velocity embedding stream plot without saving it immediately
    mv.velocity_embedding_stream(adata_result, basis='umap', color=CELL_TYPE_COLNAME)

    # Add the title
    plt.title(f"{SAMPLE_OF_INTEREST}_{CELL_TYPE_COLNAME}")

    # Save the plot to the specified file path
    plt.savefig(file_path)

    # Close the figure to free memory
    plt.close()
    
    # Confirm the plot has been saved
    print(f"Velocity_streams_per_individuals plot saved to: {file_path}")
    
    
    # Generate and save the velocity embedding stream plot
    #mv.velocity_embedding_stream(adata_result, basis='umap', color=CELL_TYPE_COLNAME, save=file_path)

    # Confirm the plot has been saved
    #print(f"Velocity stream plot saved to: {file_path}")
    
    # Construct the file path for saving the plot
    file_path = os.path.join(
        PATH_ANALYSIS_OUTPUT, 
        "Figures/Velocity_streams_per_pseudotime", 
        f"{SAMPLE_OF_INTEREST}_Velostream_pseudotime.pdf"
    )

    # Create the directories if they don't exist
    os.makedirs(os.path.dirname(file_path), exist_ok=True)

    scv.pl.scatter(adata_result, color='latent_time', color_map='gnuplot', size=80, save=file_path)

    # Confirm the plot has been saved
    print(f"Velocity stream by pseudotime plot saved to: {file_path}")
    
    
    
    #Examine some example genes downstream
    # Find the intersection of GENE_LIST_OF_INTEREST and adata_result.var_names
    gene_list = list(set(GENE_LIST_OF_INTEREST).intersection(adata_result.var_names))
    
    ## Plot accessibility and expression against gene time
    import os
    import matplotlib.pyplot as plt

    def generate_plots(adata_result, gene_list, PATH_ANALYSIS_OUTPUT, SAMPLE_OF_INTEREST, CELL_TYPE_COLNAME):
        # Define all plot configurations
        plot_configs = [
            {
                "dir_name": "Accessibility_Expression_Against_Gene_Time",
                "plot_type": "dynamic",
                "params": {
                    "color_by": "state",
                    "axis_on": False,
                    "frame_on": False
                }
            },
            {
                "dir_name": "Velocity_Against_Gene_Time",
                "plot_type": "dynamic",
                "params": {
                    "color_by": "state",
                    "by": "velocity",
                    "axis_on": False,
                    "frame_on": False
                }
            },
            {
                "dir_name": "Accessibility_Expression_Against_Shared_Latent_Time",
                "plot_type": "dynamic",
                "params": {
                    "color_by": CELL_TYPE_COLNAME,
                    "gene_time": False,
                    "axis_on": False,
                    "frame_on": False
                }
            },
            {
                "dir_name": "unspliced_unspliced",
                "plot_type": "scatter",
                "params": {
                    "color_by": CELL_TYPE_COLNAME,
                    "by": "us",
                    "axis_on": False,
                    "frame_on": False
                }
            },
            {
                "dir_name": "Unspliced_spliced",
                "plot_type": "scatter",
                "params": {
                    "color_by": "c",
                    "by": "us",
                    "cmap": "coolwarm",
                    "title_more_info": True,
                    "axis_on": False,
                    "frame_on": False
                }
            },
            {
                "dir_name": "Chromatin_unspliced",
                "plot_type": "scatter",
                "params": {
                    "color_by": CELL_TYPE_COLNAME,
                    "by": "cu",
                    "axis_on": False,
                    "frame_on": False
                }
            },
            {
                "dir_name": "3Dphaseportraits",
                "plot_type": "scatter",
                "params": {
                    "color_by": CELL_TYPE_COLNAME,
                    "by": "cus",
                    "axis_on": False,
                    "downsample": 2
                }
            },
            {
                "dir_name": "Arrow1",
                "plot_type": "scatter",
                "params": {
                    "color_by": CELL_TYPE_COLNAME,
                    "by": "us",
                    "axis_on": False,
                    "frame_on": False,
                    "downsample": 2,
                    "velocity_arrows": True
                }
            },
            {
                "dir_name": "Arrow2",
                "plot_type": "scatter",
                "params": {
                    "color_by": CELL_TYPE_COLNAME,
                    "by": "cu",
                    "axis_on": False,
                    "frame_on": False,
                    "downsample": 2,
                    "velocity_arrows": True
                }
            },
            {
                "dir_name": "Arrow3",
                "plot_type": "scatter",
                "params": {
                    "color_by": CELL_TYPE_COLNAME,
                    "by": "cus",
                    "downsample": 3,
                    "velocity_arrows": True
                }
            }
        ]

        # Generate each plot
        for config in plot_configs:
            # Construct file path
            file_path = os.path.join(
                PATH_ANALYSIS_OUTPUT,
                "Figures/Downstream_Genes_Analysis",
                config["dir_name"],
                f"{SAMPLE_OF_INTEREST}_{config['dir_name']}.pdf"
            )

            # Create directory if it doesn't exist
            os.makedirs(os.path.dirname(file_path), exist_ok=True)

            # Generate plot based on type
            if config["plot_type"] == "dynamic":
                mv.dynamic_plot(adata_result, gene_list, **config["params"])
            else:  # scatter plot
                mv.scatter_plot(adata_result, gene_list, **config["params"])

            # Save and close figure
            plt.savefig(file_path, format='pdf')
            plt.close()
        
    # Call the function with appropriate arguments    
    generate_plots(
        adata_result,
        gene_list,
        PATH_ANALYSIS_OUTPUT,
        SAMPLE_OF_INTEREST,
        CELL_TYPE_COLNAME
    )