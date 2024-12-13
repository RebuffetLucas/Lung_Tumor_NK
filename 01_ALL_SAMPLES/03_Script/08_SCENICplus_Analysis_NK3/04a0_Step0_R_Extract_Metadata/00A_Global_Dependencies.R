# ##################################################
# Global declarations and libraries for the analysis
# ##################################################


######## R Libraries
library(Matrix)
library( digest)        # digest (hashing)
library( DT)            # datatable
library( forcats)       # fct_inorder (factors)
library( fs)            # path_sanitize
library( future)        # plan (multicore)
library( ggplot2)
#library( ggpubr)
library( ggrepel)
#library( grid)
#library( gridExtra)
library( htmltools)     # browsable
library( htmlwidgets)   # JS (for datatable)
#library( iheatmapr)     # iheatmap
library( knitr)
#library( magrittr)
library( pander)        # pander
library( patchwork)     # +/ (ggplot layout)
library( pheatmap)      # pheatmap
library( plotly)
library( plyr)
library( dplyr)
library( rmarkdown)
library( scales)        # hue_pal

# Single-cell technology
library( Seurat)

#Bertrand
library(gghalves)
library(ggdist)
library(patchwork)
library(ComplexHeatmap)


#library(reshape2)
#library(rlist)

#Load library
#library(cisTopic)
library(Seurat)
library(SeuratObject)
#library(SeuratDisk)
library(arrow)
library(Signac)
