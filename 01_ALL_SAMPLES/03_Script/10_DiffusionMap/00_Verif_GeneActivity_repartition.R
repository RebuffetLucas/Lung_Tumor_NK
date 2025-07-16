# Load required libraries
library(ggplot2)
library(reshape2)
library(dplyr)
library(patchwork)

# Load objects
PBMC <- readRDS(file.path(PATH_EXPERIMENT_RAWDATA, "RDS_FILES", SAMPLE_NK3_Only_Preprocessed))

set.seed(SEED)

# RNA Data
DefaultAssay(PBMC) <- "RNA"
var.features <- VariableFeatures(PBMC)
count.data <- FetchData(object = PBMC, assay = "RNA", layer = "counts", vars = var.features)
log.norm.data <- FetchData(object = PBMC, assay = "RNA", layer = "data", vars = var.features)

# Gene Activity Data
DefaultAssay(PBMC) <- "GeneActivity"
PBMC <- FindTopFeatures(PBMC, min.cutoff = "q95")
var.features.gene_activity <- VariableFeatures(PBMC)
count.data.gene_activity <- FetchData(object = PBMC, assay = "GeneActivity", layer = "counts", vars = var.features.gene_activity)
log.norm.data.gene_activity <- FetchData(object = PBMC, assay = "GeneActivity", layer = "data", vars = var.features.gene_activity)

########################################
#### Working per unique cells  #########
########################################

# Convert matrices to data frames
df_gene_activity <- as.data.frame(log.norm.data.gene_activity)
df_rna <- as.data.frame(log.norm.data)

# Convert rownames (cell barcodes) into an explicit column
df_gene_activity$Cell <- rownames(df_gene_activity)
df_rna$Cell <- rownames(df_rna)

# Melt data, using "Cell" as an ID variable
df_gene_activity_long <- melt(df_gene_activity, id.vars = "Cell", variable.name = "Gene", value.name = "Expression")
df_rna_long <- melt(df_rna, id.vars = "Cell", variable.name = "Gene", value.name = "Expression")

# Add source column
df_gene_activity_long$Dataset <- "Gene Activity"
df_rna_long$Dataset <- "RNA Expression"

# Combine both datasets
df_combined <- rbind(df_gene_activity_long, df_rna_long)

# **1. Violin + Box Plot (No Jitter, Optimized)**
p1 <- ggplot(df_combined, aes(x = Dataset, y = Expression, fill = Dataset)) +
  geom_violin(alpha = 0.6, color = "black", trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "black") +
  theme_minimal() +
  labs(title = "Comparison of Gene Activity vs RNA Expression",
       y = "Expression Value", x = "") +
  scale_fill_manual(values = c("Gene Activity" = "blue", "RNA Expression" = "red")) +
  theme(legend.position = "none")

# **2. Density Plot (Optimized for Large Data)**
p2 <- ggplot(df_combined, aes(x = Expression, fill = Dataset)) +
  geom_density(alpha = 0.5, adjust = 1.5) +
  theme_minimal() +
  labs(title = "Density Plot of Gene Activity vs RNA Expression",
       x = "Expression Value", y = "Density") +
  scale_fill_manual(values = c("Gene Activity" = "blue", "RNA Expression" = "red")) +
  theme(legend.position = "top")

# **3. Log-Scaled Density Plot (Avoiding Zero Overload)**
p3 <- ggplot(df_combined %>% filter(Expression > 0), aes(x = log1p(Expression), fill = Dataset)) +
  geom_density(alpha = 0.5, adjust = 1.5) +
  theme_minimal() +
  labs(title = "Log-scaled Density Plot (Expression > 0)",
       x = "Log(Expression + 1)", y = "Density") +
  scale_fill_manual(values = c("Gene Activity" = "blue", "RNA Expression" = "red")) +
  theme(legend.position = "top")

# **Save Plots**
#dir.create(file.path(PATH_ANALYSIS_OUTPUT, "Verif_GeneActivity_Slot"), showWarnings = FALSE) # Ensure directory exists
pdf(file.path(PATH_ANALYSIS_OUTPUT, "Verif_GeneActivity_Slot", "DiagnosticPlots.pdf"),
    width = 24 / 2.54, height = 24 / 2.54)
print(p1 + p2 + p3 + plot_layout(ncol = 2))
dev.off()


# Compute summary statistics
summary_stats <- df_combined %>%
  group_by(Dataset) %>%
  summarise(
    Mean = mean(Expression, na.rm = TRUE),
    Median = median(Expression, na.rm = TRUE),
    Min = min(Expression, na.rm = TRUE),
    Max = max(Expression, na.rm = TRUE),
    SD = sd(Expression, na.rm = TRUE),
    .groups = "drop"
  )

# Print summary statistics
print(summary_stats)



########################################
#### Working per genes together  #########
########################################

# Load required libraries
library(ggplot2)
library(reshape2)
library(dplyr)
library(patchwork)

# Convert matrices to data frames
df_gene_activity <- as.data.frame(log.norm.data.gene_activity)
df_rna <- as.data.frame(log.norm.data)

# Compute median expression per gene
median_gene_activity <- data.frame(
  Gene = colnames(df_gene_activity),
  Expression = apply(df_gene_activity, 2, median, na.rm = TRUE),
  Dataset = "Gene Activity"
)

median_rna <- data.frame(
  Gene = colnames(df_rna),
  Expression = apply(df_rna, 2, median, na.rm = TRUE),
  Dataset = "RNA Expression"
)

# Combine datasets
df_combined <- rbind(median_gene_activity, median_rna)

# **1. Box Plot of Median Expression per Gene**
p1 <- ggplot(df_combined, aes(x = Dataset, y = Expression, fill = Dataset)) +
  geom_violin(alpha = 0.6, color = "black", trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.2, outlier.shape = NA, color = "black") +
  theme_minimal() +
  labs(title = "Comparison of Median Gene Activity vs RNA Expression",
       y = "Median Expression Value", x = "") +
  scale_fill_manual(values = c("Gene Activity" = "blue", "RNA Expression" = "red")) +
  theme(legend.position = "none")

# **2. Density Plot of Median Expression per Gene**
p2 <- ggplot(df_combined, aes(x = Expression, fill = Dataset)) +
  geom_density(alpha = 0.5, adjust = 1.5) +
  theme_minimal() +
  labs(title = "Density Plot of Median Gene Activity vs RNA Expression",
       x = "Median Expression Value", y = "Density") +
  scale_fill_manual(values = c("Gene Activity" = "blue", "RNA Expression" = "red")) +
  theme(legend.position = "top")

# **3. Scatter Plot of Gene Activity vs RNA Expression Medians**
p3 <- ggplot(df_combined %>% spread(Dataset, Expression), aes(x = `Gene Activity`, y = `RNA Expression`)) +
  geom_point(alpha = 0.5, color = "blue") +
  geom_smooth(method = "lm", color = "red", linetype = "dashed") +  # Linear regression trendline
  theme_minimal() +
  labs(title = "Scatter Plot of Median Gene Activity vs RNA Expression",
       x = "Median Gene Activity", y = "Median RNA Expression")

# **4. Log-Scaled Scatter Plot**
p4 <- ggplot(df_combined %>% spread(Dataset, Expression), aes(x = log1p(`Gene Activity`), y = log1p(`RNA Expression`))) +
  geom_point(alpha = 0.5, color = "blue") +
  geom_smooth(method = "lm", color = "red", linetype = "dashed") +
  theme_minimal() +
  labs(title = "Log-Scaled Scatter Plot of Median Gene Activity vs RNA Expression",
       x = "Log(Median Gene Activity + 1)", y = "Log(Median RNA Expression + 1)")

# **Save Plots**
pdf(file.path(PATH_ANALYSIS_OUTPUT, "Verif_GeneActivity_Slot", "PerGene_Median_DiagnosticPlots.pdf"),
    width = 24 / 2.54, height = 24 / 2.54)
print(p1 + p2 + p3 + p4 + plot_layout(ncol = 2))
dev.off()
