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
#library( pander)        # pander
library( patchwork)     # +/ (ggplot layout)
library( pheatmap)      # pheatmap
library( plotly)
#library( plyr)
#library( dplyr)
library( rmarkdown)
library( scales)        # hue_pal

# Single-cell technology
library( Seurat)

# Functional Enrichment analysis
#library( biomaRt)
#library( clusterProfiler)
#library( org.Mm.eg.db)
#library( org.Hs.eg.db)
#library( rrvgo)

#Bertrand
#library(gghalves)
#library(ggdist)
library(patchwork)
library(ComplexHeatmap)

#Me
library(dplyr)
library(plyr)
library(readxl)
library(RColorBrewer)


#For batch correction with Seurat
library(BiocManager)
#library(multtest)
#library(metap)

#Harmony
library(harmony)
library(batchelor)
library(limma)
#library(unix)
library(DESeq2)
library(rlist)
library(reshape2)
library(rlist)

library(SeuratDisk)

#Def function
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))}



#Install Libraries
#Works
#install.packages("sinaplot")
#BiocManager::install("destiny")

#install.packages("rgl")
#install.packages("magick")
#BiocManager::install("genefilter")

#Does not work
#remotes::install_github("Japrin/sscVis")
#devtools::install_github("Japrin/sscVis")



#Library Import

library(Seurat)
library(dplyr)
#library(sscVis)
library(harmony)
#library(UCell)
library(patchwork)
#library(Matrix.utils)
library(DESeq2)
library(ggplot2)
library(RColorBrewer)
#library(sinaplot)
library(destiny)
#library(rgl)
#library(magick)
library(genefilter)
#library(SeuratWrappers)
library("scatterplot3d") 
library(batchelor)


#remotes::install_github('satijalab/seurat-wrappers')

# Functions
exp2col<-function(value,breaks=1000,
                  cols=viridis::plasma(11, alpha=1, begin = 0, end=1, direction = 1)){
  colPal <- colorRampPalette(cols)
  marker.color <- colPal(breaks)[as.numeric(cut(value, breaks = breaks))]
  return(marker.color)
}



