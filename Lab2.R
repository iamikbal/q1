# Question 1:
protein <- read.csv("eur_protein_consump.csv", row.names = 1)
protein <- scale(protein)
# (a) Obtain principal component projection of the observed data.
# Keep the case tag on the projected plane.

# Perform PCA
pca_protein <- prcomp(protein,
                      center = TRUE,
                      scale. = TRUE)

summary(pca_protein)
# Extract Principal Component Scores
scores <- pca_protein$x
head(scores)

# 2D PCA Projection (PC1 vs PC2 with country labels)
plot(scores[,1], scores[,2],
     xlab = "Principal Component 1",
     ylab = "Principal Component 2",
     main = "PCA Projection (EU Protein Data)",
     pch = 16)

text(scores[,1], scores[,2],
     labels = rownames(scores),
     pos = 3, cex = 0.7)


# (b) Proportion of total variation explained
# + cumulative proportion plot

# Eigenvalues
eigenvalues <- (pca_protein$sdev)^2

# Proportion of variance explained
prop_var <- eigenvalues / sum(eigenvalues)
prop_var
# Cumulative proportion
cum_var <- cumsum(prop_var)
cum_var

# Cumulative Variance Plot
plot(cum_var,
     type = "b",
     pch = 19,
     xlab = "Number of Principal Components",
     ylab = "Cumulative Proportion of Variance",
     main = "Cumulative Explained Variance")

abline(h = 0.9, col = "red", lty = 2)

# (c) Scree Plot and Elbow Method
plot(eigenvalues,
     type = "b",
     pch = 16,
     xlab = "Principal Component",
     ylab = "Eigenvalue",
     main = "Scree Plot")
abline(h = 1, col="red", lty = 2)
# (d) Identify Multidimensional Outliers
# Use Mahalanobis distance on first 2 PCs.
pc12 <- scores[,1:2]

md <- mahalanobis(pc12,
                  colMeans(pc12),
                  cov(pc12))

cutoff <- qchisq(0.975, df = 2)

outliers <- which(md > cutoff)

outliers
rownames(scores)[outliers]

# Highlight Outliers on Plot
plot(pc12[,1], pc12[,2],
     xlab = "PC1",
     ylab = "PC2",
     main = "Outlier Detection in PCA Space",
     pch = 19)

text(pc12[,1], pc12[,2],
     labels = rownames(scores),
     pos = 3, cex = 0.7)

points(pc12[outliers,1],
       pc12[outliers,2],
       col = "red",
       pch = 19)

# (e) Identify Clusters and Their Characteristics
# Perform k-means clustering in PCA space.
set.seed(123)
kmeans_result <- kmeans(pc12, centers = 3)
clusters <- kmeans_result$cluster

# Plot Clusters
plot(pc12[,1], pc12[,2],
     col = clusters,
     pch = 19,
     xlab = "PC1",
     ylab = "PC2",
     main = "Clusters in PCA Space")

text(pc12[,1], pc12[,2],
     labels = rownames(scores),
     pos = 3, cex = 0.7)
# Cluster Characteristics (Mean Protein Intake)
aggregate(protein,
          by = list(Cluster = clusters),
          mean)

# Question 2: 
bank_raw <- read.csv("PS_bank_fin_ratio.csv",
                     stringsAsFactors = FALSE)

# Save bank names (first column, excluding header row)
bank_names <- bank_raw[-1, 1]

# Keep only numeric part (remove first row and first column)
bank_clean <- as.data.frame(lapply(bank_raw[-1, -1], as.numeric))

# Assign bank names as row names
rownames(bank_clean) <- bank_names

# Check structure
str(bank_clean)


# (a) Obtain principal component projection of the observed data.
# Keep case tag on projected plane.
# Perform PCA
pca_bank <- prcomp(bank_clean,
                   center = TRUE,
                   scale. = TRUE)

summary(pca_bank)

# Extract PC Scores
scores_bank <- pca_bank$x
# 2D PCA Projection (PC1 vs PC2 with Bank Names)

plot(scores_bank[,1], scores_bank[,2],
     xlab = "Principal Component 1",
     ylab = "Principal Component 2",
     main = "PCA Projection (Bank Financial Ratios)",
     pch = 19)

text(scores_bank[,1], scores_bank[,2],
     labels = rownames(scores_bank),
     pos = 3, cex = 0.7)

# (b) Proportion of Variance Explained
# Eigenvalues
eigenvalues_bank <- (pca_bank$sdev)^2
eigenvalues_bank
# Proportion
prop_var_bank <- eigenvalues_bank / sum(eigenvalues_bank)
prop_var_bank
# Cumulative Proportion
cum_var_bank <- cumsum(prop_var_bank)
cum_var_bank
# Cumulative Plot
plot(cum_var_bank,
     type = "b",
     pch = 19,
     xlab = "Number of Principal Components",
     ylab = "Cumulative Proportion of Variance",
     main = "Cumulative Explained Variance")

abline(h = 0.9, col = "red", lty = 2)

# (c) Scree Plot & Number of PCs to Retain
plot(eigenvalues_bank,
     type = "b",
     pch = 19,
     xlab = "Principal Component",
     ylab = "Eigenvalue",
     main = "Scree Plot")
which(eigenvalues_bank > 1)

# (d) Multidimensional Outliers
pc12_bank <- scores_bank[,1:2]

md_bank <- mahalanobis(pc12_bank,
                       colMeans(pc12_bank),
                       cov(pc12_bank))

cutoff_bank <- qchisq(0.975, df = 2)

outliers_bank <- which(md_bank > cutoff_bank)

rownames(scores_bank)[outliers_bank]

# Highlight Outliers
plot(pc12_bank[,1], pc12_bank[,2],
     xlab = "PC1",
     ylab = "PC2",
     main = "Outlier Detection",
     pch = 19)

text(pc12_bank[,1], pc12_bank[,2],
     labels = rownames(scores_bank),
     pos = 3, cex = 0.6)

points(pc12_bank[outliers_bank,1],
       pc12_bank[outliers_bank,2],
       col = "red",
       pch = 19)
# (e) Clustering and Interpretation
# K-means (choose k = 3 initially)
set.seed(123)
kmeans_bank <- kmeans(pc12_bank, centers = 3)

clusters_bank <- kmeans_bank$cluster
# Cluster Plot
plot(pc12_bank[,1], pc12_bank[,2],
     col = clusters_bank,
     pch = 19,
     xlab = "PC1",
     ylab = "PC2",
     main = "Clusters in PCA Space")

text(pc12_bank[,1], pc12_bank[,2],
     labels = rownames(scores_bank),
     pos = 3, cex = 0.6)
# Cluster Characteristics
aggregate(bank_clean,
          by = list(Cluster = clusters_bank),
          mean)

