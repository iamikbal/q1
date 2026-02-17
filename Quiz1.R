# Read data
data <- read.csv("cyber_crime.csv", header = TRUE)

# View structure
str(data)

# Set State/UT names as row names (assuming first column is State/UT)
rownames(data) <- data[,1]

# Keep only numeric columns
data_numeric <- data[,-1]
data <- data[, sum(is.na(data)) > 0]

# Check missing values
sum(is.na(data_numeric))
is.na(data_numeric)
# Standardize data
data_scaled <- scale(data_numeric)


# (a) PCA Projection
# Perform PCA
pca_result <- prcomp(data_scaled, center = TRUE, scale. = TRUE)

# Summary
summary(pca_result)

# PCA Scores (projection)
pca_scores <- pca_result$x
head(pca_scores)

# Scatter Plot of First Two PCs
plot(pca_scores[,1], pca_scores[,2],
     xlab = "PC1",
     ylab = "PC2",
     main = "PCA Projection",
     pch = 19)

text(pca_scores[,1], pca_scores[,2],
     labels = rownames(data_scaled),
     pos = 3, cex = 0.6)
pc12 <- pca_scores[,1:16]
md <- mahalanobis(pc12, colMeans(pc12), cov(pc12))
cutoff <- qchisq(0.975, df = ncol(data_scaled))
outliers <- which(md > cutoff)
rownames(data_scaled)[outliers]

plot(pc12[,1], pc12[,2],
     xlab = "PC1",
     ylab = "PC2",
     main = "Outlier detection in PC space",
     pch = 19)
text(pc12[,1], pc12[,2],
     labels = rownames(data_scaled),
     pos = 3, cex = 0.6)
points(pc12[outliers, 1], pc12[outliers, 2],
       col = "red",
       pch = 19)
# (b) Outlier Detection (Mahalanobis Distance)
# Compute Mahalanobis distance
center <- colMeans(data_scaled)
cov_matrix <- cov(data_scaled)

md <- mahalanobis(data_scaled, center, cov_matrix)

# Chi-square cutoff (97.5%)
cutoff <- qchisq(0.975, df = ncol(data_scaled))

# Identify outliers
outliers <- which(md > cutoff)

# Print outlier states
rownames(data_scaled)[outliers]

plot(md, type = "h",
     main = "Mahalanobis Distance",
     ylab = "Distance")

abline(h = cutoff, col = "red", lty = 2)

# (c) Scree Plot
# Eigenvalues
eigenvalues <- pca_result$sdev^2

# Scree Plot
plot(eigenvalues, type = "b",
     xlab = "Principal Component",
     ylab = "Eigenvalue",
     main = "Scree Plot")

abline(h = 1, col = "red", lty = 2)

# (d) Complete Linkage Hierarchical Clustering
# Distance matrix
dist_matrix <- dist(data_scaled)

# Complete linkage clustering
hc_complete <- hclust(dist_matrix, method = "complete")

# Plot dendrogram
plot(hc_complete,
     main = "Complete Linkage Dendrogram",
     cex = 0.6)

# Height where tree is cut for 4 clusters
h_cut <- sort(hc_complete$height, decreasing = TRUE)[4]

# Add horizontal cutting line
abline(h = h_cut, col = "red", lty = 2, lwd = 2)

# Get 4 clusters
clusters_hc <- cutree(hc_complete, k = 4)

# List states in each cluster
split(rownames(data_scaled), clusters_hc)


# (e) K-Means Clustering (k = 4)
set.seed(123)

kmeans_result <- kmeans(data_scaled,
                        centers = 4,
                        nstart = 25)

# Cluster assignment
kmeans_result$cluster

# List states per cluster
split(rownames(data_scaled),
      kmeans_result$cluster)

# Plot K-means Clusters (Using First 2 PCs)
plot(pca_scores[,1], pca_scores[,2],
     col = kmeans_result$cluster,
     pch = 19,
     xlab = "PC1",
     ylab = "PC2",
     main = "K-means Clustering (k=4)")

legend("bottomleft",
       legend = 1:4,
       col = 1:4,
       pch = 19)

