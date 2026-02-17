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

