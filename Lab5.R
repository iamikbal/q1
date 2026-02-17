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
