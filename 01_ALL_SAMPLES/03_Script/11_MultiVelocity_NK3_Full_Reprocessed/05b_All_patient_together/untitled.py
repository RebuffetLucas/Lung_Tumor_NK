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

for SAMPLE_OF_INTEREST in SAMPLES_NAMES_LIST :
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
        
        
    generate_plots(
        adata_result,
        gene_list,
        PATH_ANALYSIS_OUTPUT,
        SAMPLE_OF_INTEREST,
        CELL_TYPE_COLNAME
    )