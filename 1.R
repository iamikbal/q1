###### Lab1:

# install.packages(c("GGally", "corrplot", 
#                   "reshape2", "aplpack", "gridExtra"))

library(ggplot2)
library(GGally)
library(corrplot)
library(reshape2)
library(aplpack)
library(gridExtra)

## Problem 1:
# Load the data:

eco <- read.csv("eco_dev_data.csv")

eco <- eco[-1, ]
eco_numeric <- eco[, -1]   # remove country column

eco_numeric <- as.data.frame(
  lapply(eco_numeric, function(x) as.numeric(as.character(x)))
)

colnames(eco_numeric) <- c("GNPPER", "GDPGR", "DOMINV",
                           "GDPDFL", "AGRVLAD", "INDVLAD",
                           "EXP", "GOVCON", "RESBL",
                           "DOMCRDT", "GRIIMP",
                           "IMPCOV", "INTSPRD")


# (i) Box plot:

boxplot(eco_numeric,
        las=2,
        pch = 16,
        col="lightblue",
        main="Boxplots of Economic Indicators")

# (i) IQR values:

IQR_values <- apply(eco_numeric, 2, IQR)

# max IQR:

max_IQR_var <- names(which.max(IQR_values))
max_IQR_value <- max(IQR_values)

max_iqr <- c(max_IQR_var, max_IQR_value)

# min IQR:

min_IQR_var <- names(which.min(IQR_values))
min_IQR_value <- min(IQR_values)

min_iqr <- c(min_IQR_var, min_IQR_value)

# (ii) Scatterplot matrix:

pairs(eco_numeric,
      col = rgb(0, 0, 1, 0.4),  # semi-transparent blue
      pch = 19,
      cex = 0.6,
      main = "Scatterplot Matrix of Economic Indicators")

## Outlier detection:
outliers_list <- lapply(eco_numeric, function(x) boxplot.stats(x)$out)

outliers_list

## Multivariate outlier detection:
center <- colMeans(eco_numeric)
cov_matrix <- cov(eco_numeric)
md <- mahalanobis(eco_numeric, center, cov_matrix)
p <- ncol(eco_numeric)
cutoff <- qchisq(0.975, df = p)
outlier_index <- which(md > cutoff)
eco[outlier_index, 1]

plot(md,
     pch = 19,
     col = ifelse(md > cutoff, "red", "blue"),
     main = "Mahalanobis Distance",
     xlab = "Countries",
     ylab = "Distance")
abline(h = cutoff, col = "darkgreen", lwd = 2)

## (iii) Mean Vector:
mean_vec <- colMeans(eco_numeric)

## (iv) Profile plot:
eco_scaled <- scale(eco_numeric)

par(mar = c(9,4,4,2))

matplot(t(eco_scaled),
        type = "l",
        lty = 1,
        col = rgb(0, 0, 1, 0.2),
        xaxt = "n",
        xlab = "Economic Indicators",
        ylab = "Standardized Values",
        main = "Profile Plot of Countries")

axis(1,
     at = 1:ncol(eco_scaled),
     labels = colnames(eco_scaled),
     las = 2,
     cex.axis = 0.7)

abline(h = 0, col = "red", lwd = 2)

## (v) Developed and Undeveloped Profile plots and mean vectors:
eco_scaled <- scale(eco_numeric)
developed_scaled <- eco_scaled[eco_numeric$GNPPER > 1.5, ]
underdeveloped_scaled <- eco_scaled[eco_numeric$GNPPER < -0.8, ]
mean_dev <- colMeans(developed_scaled)
mean_undev <- colMeans(underdeveloped_scaled)

group_means <- rbind(mean_dev, mean_undev)

dim(group_means)


x_vals <- 1:ncol(group_means)

par(mar = c(9,4,4,2))

matplot(x_vals,
        t(group_means),
        type = "l",
        lty = 1,
        lwd = 2,
        col = c("blue", "red"),
        xaxt = "n",
        xlab = "Economic Indicators",
        ylab = "Standardized Mean Values",
        main = "Profile Plot: Developed vs Underdeveloped")

axis(1,
     at = x_vals,
     labels = colnames(eco_scaled),
     las = 2,
     cex.axis = 0.7)

abline(h = 0, col = "darkgreen", lwd = 2)

legend("top",
       legend = c("Dev", "Undev"),
       col = c("blue", "red"),
       lwd = 2)

## (vi) Covariance and Correlation Matrix:
cov_matrix <- cov(eco_numeric)
cov_matrix <- round(cov_matrix, 3)
corr_matrix <- cor(eco_numeric)
corr_matrix <- round(corr_matrix, 3)

## (vii) Correlation Heatmap:
corrplot(corr_matrix,
         method = "color",
         type = "upper",
         tl.cex = 0.7,
         tl.col = "black",
         addCoef.col = "black",
         number.cex = 0.5)

weak_count <- rowSums(abs(corr_matrix) < 0.2) - 1  # subtract 1 for diagonal
names(weak_count[weak_count > (ncol(corr_matrix)/2)])

## (viii) Maxm and Minm Correlation: 
corr_no_diag <- corr_matrix
diag(corr_no_diag) <- NA

max_value <- max(corr_no_diag, na.rm = TRUE)

max_index <- which(corr_no_diag == max_value, arr.ind = TRUE)

max_value
max_index
var1_max <- colnames(corr_no_diag)[max_index[1]]
var2_max <- colnames(corr_no_diag)[max_index[2]]

var1_max
var2_max

min_value <- min(corr_no_diag, na.rm = TRUE)

min_index <- which(corr_no_diag == min_value, arr.ind = TRUE)

min_value
min_index

var1_min <- colnames(corr_no_diag)[min_index[1]]
var2_min <- colnames(corr_no_diag)[min_index[2]]

var1_min
var2_min

## (ix) Chernoff faces:
eco_scaled <- scale(eco_numeric)
faces(eco_scaled)


##### Lab 2:

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

#### Lab3:

# ==========================================
# MTH 443 – Lab Problem Set 3
# Dataset 1: eco_dev_data.csv
# ==========================================


# ------------------------------------------
# 1. Read Data
# ------------------------------------------
data_raw <- read.csv("eco_dev_data.csv", stringsAsFactors = FALSE)

# First column assumed to be country names
# Convert all other columns to numeric
data_raw[,-1] <- lapply(data_raw[,-1], function(x) as.numeric(as.character(x)))

# Remove rows with missing values
data_clean <- na.omit(data_raw)

# Separate country names and numeric data
country_names <- data_clean[,1]
data <- data_clean[,-1]
rownames(data) <- country_names
# Check dimensions
dim(data)

# ------------------------------------------
# 2. Standardize Data
# ------------------------------------------
data_scaled <- scale(data)

# ------------------------------------------
# 3. Perform PCA
# ------------------------------------------
pca_result <- prcomp(data_scaled, center = TRUE, scale. = TRUE)

# PCA Summary
summary(pca_result)

# Eigenvalues
eigenvalues <- pca_result$sdev^2
eigenvalues

# Scree Plot
plot(eigenvalues, type = "b", pch = 19,
     xlab = "Principal Component",
     ylab = "Eigenvalue",
     main = "Scree Plot")

# ------------------------------------------
# (i) Feature Clustering using PCA
# ------------------------------------------

# Extract loadings
loadings <- pca_result$rotation
print(loadings)

# Biplot
biplot(pca_result, cex = 0.6)

# Hierarchical clustering of variables
# Use first 3 PCs for clustering
loadings_subset <- loadings[,1:3]

dist_matrix <- dist(loadings_subset)
hc <- hclust(dist_matrix, method = "ward.D2")

plot(hc, main = "Feature Clustering using PCA Loadings")
h_cut <- sort(hc$height, decreasing = TRUE)[3]
abline(h = h_cut, col = "red", lwd = 2, lty = 2)
clusters_hc <- cutree(hc, k = 3)
split(colnames(data), clusters_hc)
# ------------------------------------------
# (ii) PCA-Based Ranking of Countries
# ------------------------------------------

scores <- pca_result$x
ranking_index <- scores[,1]   # PC1 as development index

# Automatically detect GNPPER column
gnp_col <- which(grepl("GNPPER", colnames(data)))

# Ensure higher GNPPER => higher rank
if(length(gnp_col) > 0){
  if(loadings[gnp_col, 1] < 0){
    ranking_index <- -ranking_index
  }
}

ranking_table <- data.frame(
  Country = country_names,
  PCA_Index = ranking_index
)

# Sort descending
ranking_table <- ranking_table[order(-ranking_table$PCA_Index), ]
ranking_table$Rank <- 1:nrow(ranking_table)

# Top 10 countries
head(ranking_table, 10)

# ------------------------------------------
# (iii) Multivariate Normality Check
# Without using MVN package
# ------------------------------------------

# Use PCA scores
scores <- pca_result$x

n <- nrow(scores)
p <- ncol(scores)

# ------------------------------------------
# 1. Mahalanobis Distances
# ------------------------------------------

center <- colMeans(scores)
cov_mat <- cov(scores)

D2 <- mahalanobis(scores, center, cov_mat)

# ------------------------------------------
# 2. Chi-square Q-Q Plot
# ------------------------------------------

D2_sorted <- sort(D2)
chi_sq_quantiles <- qchisq(ppoints(n), df = p)

plot(chi_sq_quantiles, D2_sorted,
     main = "Chi-square Q-Q Plot",
     xlab = "Theoretical Quantiles",
     ylab = "Ordered Mahalanobis Distances")

abline(0, 1, col = "red", lwd = 2)

# If points lie approximately on line → MVN plausible


# ==========================================
# MTH 443 – Lab 3
# Dataset 2: PS_bank_fin_ratio.csv
# ==========================================


# ------------------------------------------
# 1. Read and Clean Data
# ------------------------------------------

data_raw <- read.csv("PS_bank_fin_ratio.csv",
                     stringsAsFactors = FALSE)

# Check names
print(colnames(data_raw))

# Assume:
# Column 1 = Bank
# Column 2 = Year
# Remaining columns = financial ratios

# Convert ratio columns to numeric
data_raw[ , -(1:2)] <- lapply(data_raw[ , -(1:2)],
                              function(x) as.numeric(as.character(x)))

# Remove rows with missing values
data_clean <- na.omit(data_raw)

# Verify dimensions
dim(data_clean)


# ------------------------------------------
# 2. PCA-Based Ranking for Each Year
# ------------------------------------------

# ------------------------------------------
# Separate years
# ------------------------------------------

years <- sort(unique(data_clean$Year))

for (yr in years) {
  
  cat("\n====================================\n")
  cat("YEAR:", yr, "\n")
  cat("====================================\n")
  
  # Subset data for that year
  data_year <- subset(data_clean, Year == yr)
  
  bank_names <- data_year$Bank
  ratio_data <- data_year[ , -(1:2)]
  
  # Standardize
  ratio_scaled <- scale(ratio_data)
  
  # PCA
  pca_year <- prcomp(ratio_scaled)
  
  # Scores and loadings
  scores_year <- pca_year$x
  loadings_year <- pca_year$rotation
  
  # Use PC1 as performance index
  ranking_index <- scores_year[,1]
  
  # Make sure higher ROE => higher rank
  roe_col <- which(grepl("ROE", colnames(ratio_data)))
  
  if(length(roe_col) > 0){
    if(loadings_year[roe_col,1] < 0){
      ranking_index <- -ranking_index
    }
  }
  
  # Create ranking table
  ranking_table <- data.frame(
    Bank = bank_names,
    PCA_Index = ranking_index
  )
  
  ranking_table <- ranking_table[order(-ranking_table$PCA_Index), ]
  ranking_table$Rank <- 1:nrow(ranking_table)
  
  print(ranking_table)
}


# ------------------------------------------
# 3. Multivariate Normality Check
# (Using full dataset PCA scores)
# ------------------------------------------

ratio_data_all <- data_clean[ , -(1:2)]
ratio_scaled_all <- scale(ratio_data_all)

pca_all <- prcomp(ratio_scaled_all)
scores_all <- pca_all$x

n <- nrow(scores_all)
p <- ncol(scores_all)

center <- colMeans(scores_all)
cov_mat <- cov(scores_all)

D2 <- mahalanobis(scores_all, center, cov_mat)

# Q-Q plot
D2_sorted <- sort(D2)
chi_sq_quantiles <- qchisq(ppoints(n), df = p)

plot(chi_sq_quantiles, D2_sorted,
     main = "Chi-square Q-Q Plot (Banks)",
     xlab = "Theoretical Quantiles",
     ylab = "Ordered Mahalanobis Distances")

abline(0,1,col="red",lwd=2)

# ------------------------------------------
# 4. Trajectory of Banks over 4 Years
# ------------------------------------------

# PCA on entire dataset for consistent projection
scores_plot <- scores_all[,1:2]

plot(scores_plot,
     type="n",
     xlab="PC1",
     ylab="PC2",
     main="Trajectory of Banks over Years")

banks <- unique(data_clean$Bank)

for (b in banks) {
  
  bank_data <- subset(data_clean, Bank == b)
  
  indices <- which(data_clean$Bank == b)
  
  lines(scores_plot[indices,1],
        scores_plot[indices,2],
        type="b", pch=19)
}

legend("topright", legend=banks, cex=0.6)



#### Lab5:
# -------------------------------
# 1️⃣ DATA CLEANING
# -------------------------------

wine <- read.csv("Wine_data.csv", stringsAsFactors = FALSE)

# Remove fully empty columns
wine <- wine[, colSums(!is.na(wine)) > 0]

# Identify numeric columns
numeric_cols <- sapply(wine, is.numeric)

# Keep only rows with complete numeric data
wine_clean <- wine[complete.cases(wine[, numeric_cols]), ]

# Define Type and feature matrix properly
type <- as.factor(wine_clean$Type)
X <- scale(wine_clean[, numeric_cols])

cat("Observations:", nrow(X), "\n")
cat("Variables:", ncol(X), "\n")

# -------------------------------
# 2️⃣ HIERARCHICAL CLUSTERING
# -------------------------------

d <- dist(X)

# ---- Complete Linkage ----
hc_complete <- hclust(d, method = "complete")

n <- nrow(X)
h_complete <- hc_complete$height[n - 2]

plot(hc_complete,
     labels = FALSE,
     hang = -1,
     main = "Complete Linkage Dendrogram")

abline(h = h_complete, col = "red", lwd = 2, lty = 2)

cl_complete <- cutree(hc_complete, k = 3)

# ---- Average Linkage ----
hc_average <- hclust(d, method = "average")

h_average <- hc_average$height[n - 2]

plot(hc_average,
     labels = FALSE,
     hang = -1,
     main = "Average Linkage Dendrogram")

abline(h = h_average, col = "blue", lwd = 2, lty = 2)

cl_average <- cutree(hc_average, k = 3)

# -------------------------------
# 3️⃣ K-MEANS (5 Initializations)
# -------------------------------

set.seed(123)

k1 <- kmeans(X, centers = 3, nstart = 1)
k2 <- kmeans(X, centers = 3, nstart = 5)
k3 <- kmeans(X, centers = 3, nstart = 10)
k4 <- kmeans(X, centers = 3, nstart = 20)
k5 <- kmeans(X, centers = 3, nstart = 50)

within_ss <- c(k1$tot.withinss,
               k2$tot.withinss,
               k3$tot.withinss,
               k4$tot.withinss,
               k5$tot.withinss)

data.frame(
  nstart = c(1,5,10,20,50),
  tot_within_ss = within_ss
)

# Choose best clustering
best_index <- which.min(within_ss)
cl_kmeans <- list(k1,k2,k3,k4,k5)[[best_index]]$cluster

# -------------------------------
# 4️⃣ TRACE(W) AND TRACE(B)
# -------------------------------

cluster_scatter <- function(X, cluster){
  
  X <- as.matrix(X)
  p <- ncol(X)
  overall_mean <- colMeans(X)
  
  W <- matrix(0, p, p)
  B <- matrix(0, p, p)
  
  for(k in unique(cluster)){
    
    Xk <- X[cluster == k, , drop = FALSE]
    nk <- nrow(Xk)
    mean_k <- colMeans(Xk)
    
    W <- W + t(Xk - matrix(mean_k, nk, p, byrow=TRUE)) %*%
      (Xk - matrix(mean_k, nk, p, byrow=TRUE))
    
    diff <- matrix(mean_k - overall_mean, p, 1)
    B <- B + nk * diff %*% t(diff)
  }
  
  list(traceW = sum(diag(W)),
       traceB = sum(diag(B)))
}

complete_scatter <- cluster_scatter(X, cl_complete)
average_scatter  <- cluster_scatter(X, cl_average)
kmeans_scatter   <- cluster_scatter(X, cl_kmeans)

results <- rbind(
  Complete = unlist(complete_scatter),
  Average  = unlist(average_scatter),
  Kmeans   = unlist(kmeans_scatter)
)

results

# -------------------------------
# 5️⃣ VALIDATION WITH TYPE
# -------------------------------

cat("\nContingency Tables:\n")
table(type, cl_complete)
table(type, cl_average)
table(type, cl_kmeans)

accuracy <- function(true, cluster){
  tab <- table(true, cluster)
  sum(apply(tab, 2, max)) / length(true)
}

cat("\nAccuracy:\n")
c(Complete = accuracy(type, cl_complete),
  Average  = accuracy(type, cl_average),
  Kmeans   = accuracy(type, cl_kmeans))


# Question 2:
crime <- read.csv("CrimesOnWomenData.csv", stringsAsFactors = FALSE)

# Remove empty columns
crime <- crime[, colSums(!is.na(crime)) > 0]

# Keep only rows with complete numeric data
numeric_cols <- sapply(crime, is.numeric)
crime_clean <- crime[complete.cases(crime[, numeric_cols]), ]

unique(crime_clean$Year)

year_data <- subset(crime_clean, Year == 2021)

states <- year_data$State

# Remove State and Year columns
X_raw <- year_data[, !(colnames(year_data) %in% c("State","Year"))]

# Standardize
X <- scale(X_raw)

nrow(X)

d <- dist(X)
n <- nrow(X)
# (i) Single Linkage
hc_single <- hclust(d, method="single")

h_single <- hc_single$height[n - 4]   # because K=5 → n-5+1 = n-4

plot(hc_single,
     labels = states,
     main = "Single Linkage Dendrogram")

abline(h = h_single, col="red", lwd=2, lty = 2)

cl_single <- cutree(hc_single, k=5)

# (ii) Average Linkage
hc_average <- hclust(d, method="average")

h_average <- hc_average$height[n - 4]

plot(hc_average,
     labels = states,
     main = "Average Linkage Dendrogram")

abline(h = h_average, col="blue", lwd=2, lty = 2)

cl_average <- cutree(hc_average, k=5)

set.seed(123)

k1 <- kmeans(X, centers=5, nstart=1)
k2 <- kmeans(X, centers=5, nstart=5)
k3 <- kmeans(X, centers=5, nstart=10)
k4 <- kmeans(X, centers=5, nstart=20)
k5 <- kmeans(X, centers=5, nstart=50)

within_ss <- c(k1$tot.withinss,
               k2$tot.withinss,
               k3$tot.withinss,
               k4$tot.withinss,
               k5$tot.withinss)

data.frame(
  nstart = c(1,5,10,20,50),
  tot_within_ss = within_ss
)

best_index <- which.min(within_ss)
cl_kmeans <- list(k1,k2,k3,k4,k5)[[best_index]]$cluster

cluster_scatter <- function(X, cluster){
  
  X <- as.matrix(X)
  p <- ncol(X)
  overall_mean <- colMeans(X)
  
  W <- matrix(0, p, p)
  B <- matrix(0, p, p)
  
  for(k in unique(cluster)){
    
    Xk <- X[cluster == k, , drop=FALSE]
    nk <- nrow(Xk)
    mean_k <- colMeans(Xk)
    
    W <- W + t(Xk - matrix(mean_k, nk, p, byrow=TRUE)) %*%
      (Xk - matrix(mean_k, nk, p, byrow=TRUE))
    
    diff <- matrix(mean_k - overall_mean, p, 1)
    B <- B + nk * diff %*% t(diff)
  }
  
  list(traceW = sum(diag(W)),
       traceB = sum(diag(B)))
}
single_scatter  <- cluster_scatter(X, cl_single)
average_scatter <- cluster_scatter(X, cl_average)
kmeans_scatter  <- cluster_scatter(X, cl_kmeans)

results <- rbind(
  Single  = unlist(single_scatter),
  Average = unlist(average_scatter),
  Kmeans  = unlist(kmeans_scatter)
)

results



### Lab 6:
# Install once if needed
# install.packages("FNN")
par(mfrow = c(1, 1))
library(FNN)
library(MASS)

chd <- read.csv("chd.csv")

# Split groups
chd1 <- subset(chd, TenYearCHD == 1)
chd0 <- subset(chd, TenYearCHD == 0)

k <- 30

# Group 1
x1 <- chd1$totChol
n1 <- length(x1)

dist1 <- get.knn(matrix(x1, ncol=1), k=k+1)$nn.dist[,k+1]
f1 <- k / (n1 * 2 * dist1)

# Group 0
x0 <- chd0$totChol
n0 <- length(x0)

dist0 <- get.knn(matrix(x0, ncol=1), k=k+1)$nn.dist[,k+1]
f0 <- k / (n0 * 2 * dist0)

plot(sort(x1), f1[order(x1)],
     type="l", col="red", lwd=2,
     main="kNN Density (k=30) - totChol",
     xlab="totChol", ylab="Density")

lines(sort(x0), f0[order(x0)],
      col="blue", lwd=2)

legend("topright",
       legend=c("CHD=1","CHD=0"),
       col=c("red","blue"), lwd=2)

k_values <- c(10, 30, 60)

par(mfrow=c(1,3))

for(k in k_values){
  dist1 <- get.knn(matrix(x1, ncol=1), k=k+1)$nn.dist[,k+1]
  f1 <- k / (n1 * 2 * dist1)
  
  plot(sort(x1), f1[order(x1)],
       type="l",
       main=paste("CHD=1, k =",k),
       xlab="totChol", ylab="Density")
}
par(mfrow = c(1,1))
# (A)(ii) Kernel Density Estimation
kernels <- c("gaussian","rectangular","triangular","epanechnikov")

par(mfrow=c(2,2))

for(kern in kernels){
  d1 <- density(chd1$totChol, kernel=kern)
  plot(d1,
       main=paste("CHD=1 -", kern),
       col="red", lwd=2)
}
par(mfrow = c(1,1))
d1 <- density(chd1$totChol, kernel="gaussian")
d0 <- density(chd0$totChol, kernel="gaussian")

plot(d1, col="red", lwd=2,
     main="Gaussian KDE - totChol",
     xlab="totChol", ylab="Density")

lines(d0, col="blue", lwd=2)

legend("topright",
       legend=c("CHD=1","CHD=0"),
       col=c("red","blue"), lwd=2)

# (D) Bivariate KDE (Gaussian Kernel)
z <- kde2d(chd$sysBP,
           chd$diaBP,
           n=60)

persp(z$x, z$y, z$z,
      theta=30, phi=30,
      col="lightblue",
      xlab="sysBP",
      ylab="diaBP",
      main="Bivariate Gaussian KDE")
contour(z,
        xlab="sysBP",
        ylab="diaBP",
        main="Contour Plot")

#### Quiz1:
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


