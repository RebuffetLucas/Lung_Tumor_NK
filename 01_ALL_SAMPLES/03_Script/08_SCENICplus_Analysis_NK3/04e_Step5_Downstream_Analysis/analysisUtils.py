# analysisUtils.py
import numpy as np
import matplotlib.pyplot as plt
import adjustText

#Manual PCA visualization 
def plot_mm_line_pca(ax, eRegulon_gene_AUC, CELL_TYPE_COLNAME, color_dict_line):
    """
    Manual PCA plot colored by cell type, with average position labels.

    Parameters:
    - ax: matplotlib axis
    - eRegulon_gene_AUC: AnnData object containing PCA coordinates and cell metadata
    - CELL_TYPE_COLNAME: column name (string) in `obs` used for coloring and labeling
    - color_dict_line: dictionary mapping cell type names to colors
    """
    texts = []
    # Scatter plot of first two PCA components
    ax.scatter(
        eRegulon_gene_AUC.obsm["X_pca"][:, 0],
        eRegulon_gene_AUC.obsm["X_pca"][:, 1],
        color = [color_dict_line[line] for line in eRegulon_gene_AUC.obs["scRNA_counts:" + CELL_TYPE_COLNAME]]
    )
    # Label cluster centers
    for line in set(eRegulon_gene_AUC.obs["scRNA_counts:" + CELL_TYPE_COLNAME]):
        line_bc_idc = np.where(eRegulon_gene_AUC.obs["scRNA_counts:" + CELL_TYPE_COLNAME] == line)[0]
        avg_x, avg_y = eRegulon_gene_AUC.obsm["X_pca"][line_bc_idc, 0:2].mean(0)
        texts.append(
            ax.text(avg_x, avg_y, line, fontweight="bold")
        )
    adjustText.adjust_text(texts)




#Manual UMAP visualization
def plot_mm_line_umap(ax, eRegulon_gene_AUC, CELL_TYPE_COLNAME, color_dict_line):
    texts = []
    # Plot UMAP
    ax.scatter(
        eRegulon_gene_AUC.obsm["X_umap"][:, 0],
        eRegulon_gene_AUC.obsm["X_umap"][:, 1],
        color = [color_dict_line[line] for line in eRegulon_gene_AUC.obs["scRNA_counts:"+ CELL_TYPE_COLNAME]]
    )
    # Plot labels
    for line in set(eRegulon_gene_AUC.obs["scRNA_counts:"+ CELL_TYPE_COLNAME]):
        line_bc_idc = np.arange(len(eRegulon_gene_AUC.obs_names))[eRegulon_gene_AUC.obs["scRNA_counts:"+ CELL_TYPE_COLNAME] == line]
        avg_x, avg_y = eRegulon_gene_AUC.obsm["X_umap"][line_bc_idc, 0:2].mean(0)
        texts.append(
            ax.text(
                avg_x,
                avg_y,
                line,
                fontweight = "bold"
            )
        )
    adjustText.adjust_text(texts)


def plot_mm_line_umap_Improved(ax, eRegulon_gene_AUC, CELL_TYPE_COLNAME, color_dict_line, point_size=1, legend_loc='center left'):
    texts = []
    
    # Plot UMAP
    scatter = ax.scatter(
        eRegulon_gene_AUC.obsm["X_umap"][:, 0],
        eRegulon_gene_AUC.obsm["X_umap"][:, 1],
        c=[color_dict_line[line] for line in eRegulon_gene_AUC.obs["scRNA_counts:" + CELL_TYPE_COLNAME]],
        s=point_size,  # Set point size
        alpha=0.7  # Add some transparency
    )
    
    # Plot labels
    for line in set(eRegulon_gene_AUC.obs["scRNA_counts:" + CELL_TYPE_COLNAME]):
        line_bc_idc = np.where(eRegulon_gene_AUC.obs["scRNA_counts:" + CELL_TYPE_COLNAME] == line)[0]
        avg_x, avg_y = eRegulon_gene_AUC.obsm["X_umap"][line_bc_idc, 0:2].mean(0)
        texts.append(
            ax.text(
                avg_x,
                avg_y,
                line,
                fontweight="bold",
                fontsize=8  # Adjust font size if needed
            )
        )
    
    # Adjust text labels to avoid overlap
    adjustText.adjust_text(texts, arrowprops=dict(arrowstyle='->', color='black'))
    
    # Add legend
    unique_lines = list(set(eRegulon_gene_AUC.obs["scRNA_counts:" + CELL_TYPE_COLNAME]))
    legend_elements = [plt.Line2D([0], [0], marker='o', color='w', 
                                  markerfacecolor=color_dict_line[line], markersize=8, label=line)
                       for line in unique_lines]
    ax.legend(handles=legend_elements, loc=legend_loc, bbox_to_anchor=(1.05, 0.5), 
              title=CELL_TYPE_COLNAME, title_fontsize='large')
    
    # Set labels and title
    ax.set_xlabel('UMAP 1')
    ax.set_ylabel('UMAP 2')
    ax.set_title('UMAP Plot')
    
    # Adjust layout to prevent clipping of labels
    plt.tight_layout()
    
    return ax
