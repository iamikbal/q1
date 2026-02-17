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

