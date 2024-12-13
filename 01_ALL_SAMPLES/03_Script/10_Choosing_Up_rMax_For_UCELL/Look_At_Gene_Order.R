# Load the necessary library
library(Matrix)
# Calculate the counts per feature (row)
feature_counts <- rowSums(sparse_matrix > 0)

# Create a data frame with feature names and their counts
feature_df <- data.frame(
  Feature = rownames(sparse_matrix),
  Counts = feature_counts
)

# Rank the features based on counts (descending order)
feature_df <- feature_df[order(-feature_df$Counts), ]
feature_df$Rank <- seq_len(nrow(feature_df))

# Specify the list of feature names you want to locate
target_features <- c("XCL1", "GZMB", "KLRC2")  # Example features

# Subset the data frame to show the ranks of the target features
ranked_features <- feature_df[feature_df$Feature %in% target_features, ]

# Display the ranked features
print(ranked_features)
