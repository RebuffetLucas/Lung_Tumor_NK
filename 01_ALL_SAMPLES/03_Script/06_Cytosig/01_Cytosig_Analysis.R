#### Installing to run a first try prior to the docker update
#BiocManager::install("dittoSeq")

##### Functions #####
bubbleZP <- function(df, celltype_col, value_col, replace = T, lower = -2, upper = 2) {
  celltypes <- unique(df[,celltype_col])
  out <- as.data.frame(matrix(nrow = length(celltypes), ncol = 3))
  colnames(out) <- c('celltype', 'z', 'p')
  out$celltype <- celltypes
  real_means <- c()
  pvalues <- c()
  
  for (aCelltype in celltypes) {
    print(aCelltype)
    real_mean <- mean(subset(df, get(celltype_col) == aCelltype)[,value_col])
    real_means <- c(real_means, real_mean)
    #pvalues <- c(pvalues, pnorm(real_mean, mean(df[,value_col]), sd(df[,value_col]), lower.tail=FALSE))
    pvalues <- c(pvalues, t.test(subset(df, get(celltype_col) == aCelltype)[,value_col], 
                                 subset(df, get(celltype_col) != aCelltype)[,value_col],
                                 alternative = 'greater')$p.value)
  }
  out$z <- real_means
  out$p <- pvalues
  out <- out %>% 
    mutate(logP = -log10(p),
           Z = rescale(z, to = c(0, 1))) %>% 
    mutate(logP = if_else(logP > 10, 10, logP))
  
  return(out)
}

NK.seu <- readRDS( file.path( PATH_EXPERIMENT_RAWDATA , "RDS_FILES",SAMPLE_NK_TUMOR) )  # Read in the data
DefaultAssay(NK.seu) = "RNA"


## 1.2 Prepare Cytosig data
# Filter genes, only keep protein-coding genes
counts = NK.seu[["RNA"]]$counts
genes <- rownames(counts)

# Convert to log2(TPM + 1) and mean centralize counts as recommended by the author (https://github.com/data2intelligence/CytoSig/issues/2)
size_factor <- 1E6/colSums(counts) # convert count to TPM
norm_counts <- sweep(counts, 2, size_factor, "*")
norm_counts <- log2(norm_counts + 1)

background <- rowMeans(norm_counts)
meanCentral_counts <- sweep(norm_counts, 1, background, "-")

# Write table
set.seed(12345)
meanCentral_rand_counts <- meanCentral_counts[,sample(ncol(meanCentral_counts))] # shuffle cells first                                              

# Define the directory where the files will be saved
output_dir <- file.path(PATH_ANALYSIS_OUTPUT, "Tables")


# Loop to save parts of the meanCentral_rand_counts matrix into separate files
for(i in 1:ceiling(ncol(meanCentral_rand_counts) / numOfCellsPerRun)) {
  print(i)
  # Select a subset of columns for the current chunk
  temp_counts <- meanCentral_rand_counts[,(numOfCellsPerRun * (i - 1) + 1):(min(numOfCellsPerRun * i, ncol(meanCentral_rand_counts)))]
  
  # Define the output file path using the specified directory
  output_file <- file.path(output_dir, paste0('new_seu_cytosig_mat_part', i, '.txt'))
  
  # Write the subset to a file
  write.table(temp_counts, output_file, quote = FALSE, sep = '\t')
}



?????????????????
?????????????????  
## 1.2 Store Cytosig results (Z-score) to the seurat object
cytosig_Z.list <- list()
for (i in 1:4) {
  print(i)
  cytosig_Z.list[[i]] <- read.csv(paste0('new_cytosig_part',i,'.Zscore'), sep = '\t')
}

all_cells <- str_replace(rownames(NK.seu@meta.data), '-1', '.1')

cytosig_Z.all.res <- do.call(cbind, cytosig_Z.list)[,all_cells]

for (cytokine in cytokines){
  print(cytokine)
  NK.seu[[paste0('cytosig_Z_',cytokine)]] <- t(cytosig_Z.all.res[cytokine,])
}










