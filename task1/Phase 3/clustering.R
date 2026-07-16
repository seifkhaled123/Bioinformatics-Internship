source("../functions.R")

pcs <- read_plink_out("../Phase\ 2/outputs/PCA.eigenvec", header = FALSE)
colnames(pcs) <- c("FID", "IID", paste0("PC", 1:10))

library(ggplot2)
library(plotly)
library(factoextra)


pca_data_2d <- pcs[, c("PC1", "PC2")]
pca_data_3d <- pcs[, c("PC1", "PC2", "PC3")]

# ---------------------------------------------------------
# TASK 2: Determine appropriate number of clusters (Elbow Method)
# ---------------------------------------------------------
# This automatically calculates the Within-Cluster Sum of Squares (WCSS) 
# and plots the elbow curve for you to justify your choice.
elbow_plot <- fviz_nbclust(pca_data_2d, kmeans, method = "wss") +
  labs(title = "Elbow Method: Optimal Clusters for Qatari PCA")
ggsave("outputs/Elbow_Method_Plot.jpg", elbow_plot, width = 7, height = 5)

# ---------------------------------------------------------
# TASK 1: Apply Clustering (k-means)
# ---------------------------------------------------------
set.seed(42) # For reproducibility

# Cluster based on 2 PCs
kmeans_2d <- kmeans(pca_data_2d, centers = 3, nstart = 25)
pcs$Cluster_2D <- as.factor(kmeans_2d$cluster)

# Cluster based on 3 PCs
kmeans_3d <- kmeans(pca_data_3d, centers = 3, nstart = 25)
pcs$Cluster_3D <- as.factor(kmeans_3d$cluster)

# ---------------------------------------------------------
# DELIVERABLES: Generate the 2D and 3D Plots
# ---------------------------------------------------------
# 1. 2D Cluster Plot (PC1 vs PC2)
plot_2d <- ggplot(pcs, aes(x = PC1, y = PC2, color = Cluster_2D)) +
  geom_point(size = 3, alpha = 0.7) +
  theme_minimal() +
  labs(title = "Population Structure (2D): PC1 vs PC2")
ggsave("outputs/PCA_2D_Clusters.jpg", plot_2d, width = 7, height = 5)

# 2. 3D Cluster Plot (PC1 vs PC2 vs PC3)
plot_3d <- plot_ly(pcs, x = ~PC1, y = ~PC2, z = ~PC3, color = ~Cluster_3D,
               type = "scatter3d", mode = "markers") %>%
  layout(title = "Population Structure (3D): PC1 vs PC2 vs PC3")

# Plotly creates interactive HTML widgets, which are much better for 3D than static JPEGs.
htmlwidgets::saveWidget(plot_3d, "outputs/PCA_3D_Clusters.html")
