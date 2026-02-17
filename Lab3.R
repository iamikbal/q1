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

