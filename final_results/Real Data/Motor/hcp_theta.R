setwd("~/Desktop/final_results/Motor")
load("simu220.08.RData")

library(fields)
library(RColorBrewer)
b_stimu <- bstimu[[1]]

b_stimu_f <- (b_stimu + t(b_stimu)) / 2
# write.csv(b_stimu_f, "hcp_tongue.csv", row.names = T)
load("simu320.08.RData")
the_stimu <- estheta

# write.csv(the_stimu, "hcp_intercept.csv", row.names = T)

# Load required libraries
library(ggplot2)
library(cluster)
library(factoextra)
library(clue)  # For Hungarian Algorithm
library(pheatmap)
library(RColorBrewer)
library(readxl)

# ========================
# Load and Preprocess Data
# ========================

# Load 90x90 connectivity matrix (skip first row if necessary)
matrix_data <- read.csv("hcp_intercept.csv", header = FALSE, skip = 1)  
matrix_data <- as.matrix(matrix_data[1:50, 2:51])  # Ensure only 90x90 is taken
library(readxl)
library(grid)
# Load true labels from AAL_FULL_NAME.xlsx
# aal_data <- read_excel("AAL_FULL_NAME.xlsx")
aaldata <- read_excel("HCP_MMP.xlsx")
aaldata$region_label <- factor(aaldata$cortex)
# Extract true labels (corresponding to the first 90 rows)
true_labels <- aaldata$region_label[185:234]  
# Ensure row names and column names of matrix match true labels
rownames(matrix_data) <- colnames(matrix_data) <- true_labels

# ========================
# Standardize Data and Determine Optimal K
# ========================

# Standardize the matrix (important for K-means)
matrix_scaled <- scale(matrix_data)

# Determine optimal number of clusters (K) using Elbow & Silhouette Methods
fviz_nbclust(matrix_scaled, kmeans, method = "wss",k.max=20)  # Elbow Method
fviz_nbclust(matrix_scaled, kmeans, method = "silhouette")  # Silhouette Method

# ========================
# Perform K-means Clustering
# ========================
set.seed(42)
k_opt <- 9  # Choose best K from elbow/silhouette analysis
kmeans_result <- kmeans(matrix_scaled, centers = k_opt, nstart = 10)

# ========================
# Match K-means Clusters with True Labels (Hungarian Algorithm)
# ========================

# Create a confusion matrix (K-means clusters vs. true labels)
conf_matrix <- table(kmeans_result$cluster, true_labels)

# Ensure confusion matrix is valid
if (nrow(conf_matrix) == 0 | ncol(conf_matrix) == 0) {
  stop("Error: Confusion matrix is empty. Check clustering and labels.")
}

# Solve Hungarian Algorithm to match clusters to true labels
conf_matrix <- as.matrix(conf_matrix)  # Ensure numeric matrix
perm <- solve_LSAP(conf_matrix, maximum = TRUE)

# Validate `perm` output
if (is.null(perm) | length(perm) == 0) {
  stop("Error: No valid assignment found. Check confusion matrix.")
}

# Assign clusters to true labels using best match
cluster_map <- setNames(levels(factor(true_labels)), perm)
mapped_clusters <- cluster_map[kmeans_result$cluster]

# ========================
# Evaluate Clustering Performance
# ========================

# Compute accuracy
accuracy <- sum(diag(conf_matrix)) / sum(conf_matrix)

# Compute purity (max cluster match in each row)
purity <- apply(conf_matrix, 1, function(x) max(x) / sum(x))

# Print results
cat("\nClustering Accuracy:", accuracy, "\n")
cat("\nCluster Purity:", mean(purity), "\n")
print(table(mapped_clusters, true_labels))

# ========================
# Visualize Results with a Heatmap
# ========================

# Reorder the matrix by cluster assignment
# Order the matrix by cluster assignment
sorted_indices <- order(mapped_clusters)  
sorted_matrix <- matrix_data[sorted_indices, sorted_indices]
sorted_clusters <- mapped_clusters[sorted_indices]  # Sorted cluster assignments

sorted_matrix <- as.matrix(sorted_matrix)
rownames(sorted_matrix) <- sorted_matrix[, 1]  # First column as row names
rownames(sorted_matrix) <- colnames(sorted_matrix)  # First column as row names
colnames(sorted_matrix) <- NULL
# sorted_matrix <- sorted_matrix[, -1]  # Remove the first column
sorted_matrix <- as.matrix(sorted_matrix)  # Convert to matrix
# rownames(row_annotation) <- rownames(sorted_matrix)

# Create a numeric index for row names (to avoid duplicates)
# rownames(sorted_matrix) <- 1:nrow(sorted_matrix)  

# Create row annotation using clusters (without duplicating names)
row_annotation <- data.frame(Cluster = factor(sorted_clusters))  # Convert clusters to factor
# rownames(row_annotation) <- rownames(sorted_matrix)  # Use simple numeric names
# plot heatmap without redundant row names
heatmap_plot <-pheatmap(
  sorted_matrix,
  cluster_rows = FALSE, cluster_cols = FALSE,  # Disable hierarchical clustering
  # annotation_row = row_annotation,  # Show cluster groups instead of full region names
  color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(100),
  show_rownames = T,  # Hide detailed row labels,
)
grid.force()
# grid.ls(print = TRUE)
seekViewport("matrix.4-3-4-3")

diagonal_blocks <- list(
  list(rows = 15:22, cols = 15:22),   # Block 2
  list(rows = 22:26, cols = 22:26),   # Block 3
  list(rows = 26:31, cols = 26:31),    # Block 4 (Bottom-right)
  list(rows = 41:50, cols = 41:50)
  # list(rows = 76:90, cols = 76:90)
)

draw_main_diagonal_blocks <- function(blocks, n) {
  for (block in blocks) {
    row_start <- (50-max(block$rows))/n
    row_end <- (50-min(block$rows))/n
    col_start <- min(block$cols)/n 
    col_end <- max(block$cols)/n 
    
    grid.rect(
      x = unit((col_start + col_end) / 2, "npc"),  # Center on X-axis (Columns)
      y = unit((row_start + row_end) / 2, "npc"),  # Center on Y-axis (Rows)
      width = unit(col_end - col_start, "npc"),    # Block width
      height = unit(row_end - row_start, "npc"),   # Block height
      gp = gpar(col = "black", lty = "dashed", lwd = 2, fill = NA)  # Red dashed border
    )
  }
}

# **Draw Diagonal Blocks Properly**
draw_main_diagonal_blocks(diagonal_blocks, nrow(sorted_matrix))

library(knitr)
library(kableExtra)

# Example data (modify based on actual values from Figure 2)
table_data <- data.frame(
  Cluster = 1:4,  # Clusters 1 to 5
  Region = I(list(
    c("Early Visual","Ventral Stream Visual","Dorsal Stream Visual","MT+ Complex and Neighboring Visual Areas"),  # Cluster 1
    c("Posterior Cingulate"),  # Cluster 2
    c("Early Visual","Dorsal Stream Visual","Superior Parietal"),  # Cluster 3
    c("Somatosensory and Motor","Paracentral Lobular and Mid Cingulate")
  ))
)

# Print table
print(table_data, row.names = FALSE)

library(knitr)
library(kableExtra)

# Convert list elements to comma-separated text
formatted_table <- table_data
formatted_table$Region <- sapply(formatted_table$Region, function(x) paste(x, collapse = ", "))
rownames(formatted_table) <- NULL
# Print without column names
kable(formatted_table, col.names = NULL, format = "html") %>%
  kable_styling(full_width = FALSE)
library(gridExtra)
library(grid)
table_plot <- tableGrob(formatted_table,rows = NULL)

# Display table
grid.newpage()
grid.draw(table_plot)