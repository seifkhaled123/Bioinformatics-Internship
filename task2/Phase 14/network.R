library(ppcor)
library(igraph)

# ==========================================
# 1. Prepare Data
# ==========================================
data <- read.csv("../Phase 11/outputs/Qatari_metabolomics_Cleaned_Scaled.csv") 

# Drop the ID and target columns. We ONLY want metabolites.
data <- data[, !(names(data) %in% c("main_id", "Diabetes"))]

# Ensure everything left is purely numeric
metabolite_matrix <- data[, sapply(data, is.numeric)] 

# ==========================================
# 2. Compute Partial Correlation Matrix
# ==========================================
cat("Computing Partial Correlations... (this might take a minute)\n")
pcor_results <- pcor(metabolite_matrix)

partial_r_matrix <- pcor_results$estimate  # The correlation weights
p_val_matrix <- pcor_results$p.value       # The statistical significance

# FIX 1: Restore the metabolite names to the rows and columns of the matrix
colnames(partial_r_matrix) <- colnames(metabolite_matrix)
rownames(partial_r_matrix) <- colnames(metabolite_matrix)

# ==========================================
# 3. Thresholding (Using FDR / Benjamini-Hochberg)
# ==========================================
cat("Applying False Discovery Rate (FDR) correction...\n")

# Adjust the p-values using BH (FDR)
p_val_vector <- as.vector(p_val_matrix)
fdr_vector <- p.adjust(p_val_vector, method = "BH")
fdr_matrix <- matrix(fdr_vector, nrow = ncol(p_val_matrix))

# Apply the thresholds (FDR adjusted P-value < 0.05 AND absolute strength >= 0.15)
adj_matrix <- partial_r_matrix
adj_matrix[fdr_matrix > 0.05 | abs(partial_r_matrix) < 0.15] <- 0 
diag(adj_matrix) <- 0 # A metabolite cannot be connected to itself

# ==========================================
# 4. Build the Network Graph in R
# ==========================================
net <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", weighted = TRUE)

# Basic R Plot (Deliverable 1)
pdf("outputs/Rplots_Network.pdf", width=10, height=10) # Save properly to a PDF
plot(net, 
     vertex.size = 3, 
     vertex.label = NA, 
     edge.width = abs(E(net)$weight) * 5, 
     main = "Metabolite Partial Correlation Network (FDR < 0.05)")
dev.off()

# ==========================================
# 5. Export for Cytoscape (Deliverable 2 & 3)
# ==========================================
# Export Edge List
edges <- as_data_frame(net, what = "edges")
write.csv(edges, "outputs/Cytoscape_Edge_List.csv", row.names = FALSE)

# Export Node List
nodes <- data.frame(id = V(net)$name)
write.csv(nodes, "outputs/Cytoscape_Node_List.csv", row.names = FALSE)

# Export underlying matrix (Deliverable 3)
write.csv(partial_r_matrix, "outputs/Partial_Correlation_Matrix.csv")

cat("Success! FDR network generated. Check Cytoscape_Edge_List.csv!\n")