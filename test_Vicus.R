# Test code of Vicus Function

## Dependent Packages Installation
if (!requireNamespace("testthat", quietly = TRUE)){
    install.packages("testthat", repos="https://cran.r-project.org")
}
library("testthat")
if (!requireNamespace("dimRed", quietly = TRUE)){
    install.packages("dimRed", repos="https://cran.r-project.org")
}
library("dimRed")

## Vicus Function (Paste your Vicus code here)

Vicus <- function(X, alpha=0.9){
    ######################################
    # Argument Check
    ######################################
    .checkVicus(X, alpha)
    # Local spectrum matrix
    V <- ...
}

.checkVicus <- function(X, alpha){
	stopifnot(is.matrix(X))
	stopifnot(is.numeric(alpha))
	stopifnot(alpha > 0)
	stopifnot(alpha < 1)
}


## Simulation Dataset
X <- matrix(runif(100*40), nrow=100, ncol=40)
X[1:25, 1:10] <- 10 * X[1:25, 1:10]
X[26:50, 11:20] <- 20 * X[26:50, 11:20]
X[51:75, 21:30] <- 30 * X[51:75, 21:30]
labelX <- rep(1:4, each=25)

## Perform Vicus and simular methods against Simulation Dataset
J=2
out.Vicus <- svd(Vicus(X))$u[, 1:J]

dat <- new("dimRedData")
dat@data <- X
out.LLE <- embed(dat, "LLE", knn = 50, ndim=J)
out.HLLE <- embed(dat, "HLLE", knn = 15, ndim=J)

plot(out.Vicus, col=labelX, pch=16, cex=2)
plot(out.LLE@data@data, col=labelX, pch=16, cex=2)
plot(out.HLLE@data@data, col=labelX, pch=16, cex=2)


## Test Input object / type
### Test I-1: Object Names

expect_identical(names(formals(Vicus)), c("X", "alpha"))

### Test I-2: X

expect_identical(as.character(formals(ICA)$X), "")

### Test I-3: alpha

expect_identical(formals(ICA)$alpha, 0.9)


## Test Output object / type
### Test O-1: Object

expect_identical(is.matrix(out.Vicus), TRUE)

## Session Information

sessionInfo()
