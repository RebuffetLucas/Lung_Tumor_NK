# Load the necessary library
library(Matrix)

# Assuming `sparse_matrix` is your dgCMatrix object
# Example: sparse_matrix <- your_sparse_matrix
sparse_matrix = PBMC_Blood_All@assays[["RNA"]]@counts

# Calculate the counts per feature (row)
feature_counts <- rowSums(sparse_matrix > 0)

# Bin the counts into intervals of 500
bins <- cut(feature_counts, breaks = seq(0, max(feature_counts), by = 100), include.lowest = TRUE, right = FALSE)

# Create a table of frequencies for the bins
binned_counts <- table(bins)

# Convert the table to a data frame for easier plotting
binned_df <- as.data.frame(binned_counts)
colnames(binned_df) <- c("Interval", "Features")

# Plot the bar plot
barplot(
  binned_df$Features,
  names.arg = binned_df$Interval,
  xlab = "Number of Counts (Binned)",
  ylab = "Number of Features",
  main = "Feature Counts Distribution (Grouped by 500s)",
  las = 2, # Rotate x-axis labels for better readability
  col = "blue",
  cex.names = 0.7 # Reduce label size for better fit
)
