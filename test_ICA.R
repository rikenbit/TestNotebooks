# Test code of ICA Function

## Dependent Packages Installation

if (!requireNamespace("testthat", quietly = TRUE)){
    install.packages("testthat", repos="https://cran.r-project.org")
}
library("testthat")

## ICA Function (Paste your ICA code here)

ICA <- function(X, J,
    algorithm=c("JADE", "AuxICA1", "AuxICA2", "IPCA", "SIMBEC", "AMUSE",
        "SOBI", "FOBI", "ProDenICA", "RICA"),
    num.iter=30, thr=1E-10, verbose=FALSE){
    ######################################
    # Argument Check
    ######################################
    .checkICA(X, J, num.iter, thr, verbose)
    algorithm <- match.arg(algorithm)
    ######################################
    # Initialization (e.g. Whiteniing)
    ######################################
    int <- .initICA(X, J, thr)
    X <- int$X
    A <- int$A
    S <- int$S
    RecError <- int$RecError
    RelChange <- int$RelChange
    ######################################
    # Iteration
    ######################################
    iter <- 1
    while ((RelChange[iter] > thr) && (iter <= num.iter)) {
        if(algorithm == "JADE"){
        	A <- .JADE(X, A, S, J)
        }
        if(algorithm == "AuxICA1"){
            A <- .AuxICA1(X, A, S, J)
        }
        if(algorithm == "AuxICA2"){
            A <- .AuxICA2(X, A, S, J)
        }
        if(algorithm == "IPCA"){
            A <- .IPCA(X, A, S, J)
        }
        if(algorithm == "SIMBEC"){
            A <- .SIMBEC(X, A, S, J)
        }
        if(algorithm == "AMUSE"){
            A <- .AMUSE(X, A, S, J)
        }
        if(algorithm == "SOBI"){
            A <- .SOBI(X, A, S, J)
        }
        if(algorithm == "FOBI"){
            A <- .FOBI(X, A, S, J)
        }
        if(algorithm == "ProDenICA"){
            A <- .ProDenICA(X, A, S, J)
        }
        if(algorithm == "RICA"){
            A <- .RICA(X, A, S, J)
        }
        # After Update
        S <- X %*% ginv(A)
        X_bar <- A %*% S
        iter <- iter + 1
        RecError[iter] <- .recError(X, X_bar)
        RelChange[iter] <- .relChange(iter, RecError)
        # Verbose
        if(verbose){
             cat(paste0(iter, " / ", num.iter,
                " |Previous Error - Error| / Error = ",
                RelChange[iter], "\n"))
        }
    }
    # Output
	names(RecError) <- c("offset", 1:(iter - 1))
    names(RelChange) <- c("offset", 1:(iter - 1))
	# Output
    list(A=A, S=S, J=J, algorithm=algorithm, num.iter=num.iter,
    thr=thr, verbose=verbose, RecError=RecError, RelChange=RelChange)
}

.JADE <- function(X, A, S, J){
    ...
}

.AuxICA1 <- function(X, A, S, J){
    ...
}

.AuxICA2 <- function(X, A, S, J){
    ...
}

.IPCA <- function(X, A, S, J){
    ...
}

.SIMBEC <- function(X, A, S, J){
    ...
}

.AMUSE <- function(X, A, S, J){
    ...
}

.SOBI <- function(X, A, S, J){
    ...
}

.FOBI <- function(X, A, S, J){
    ...
}

.ProDenICA <- function(X, A, S, J){
    ...
}

.RICA <- function(X, A, S, J){
    ...
}

.checkICA <- function(X, J, num.iter, thr, verbose){
	stopifnot(is.matrix(X))
	stopifnot(is.numeric(J))
	stopifnot(length(J) == 1)
	stopifnot(min(dim(X)) >= J)
	stopifnot(is.numeric(num.iter))
	stopifnot(num.iter > 0)
	stopifnot(is.numeric(thr))
	stopifnot(is.logical(verbose))
}

.initICA <- function(X, J, thr){
	X <- .whitening(X)
	# A/S
	nr <- nrow(X)
	A <- matrix(runif(nr*J), nrow=nr, ncol=J)
	S <- X %*% ginv(A)
	X_bar <- A %*% S
	# Reconstruction Error
	RecError <- .recError(X, X_bar)
	# Relative Change
	RelChange <- thr * 10
	list(X=X, A=A, S=S, RecError=RecError, RelChange=RelChange)
}

.whitening <- function(X){
	...
}

.recError <- function(X, X_bar){
	...
}

.relChange <- function(iter, RecError){
	...
}


## Simulation Datasets

X <- matrix(runif(100*200), nrow=100, ncol=200)

## Perform ICA against Simulation Datasets

J <- 5
out.JADE <- ICA(X, J=J, algorithm="JADE")
out.AuxICA1 <- ICA(X, J=J, algorithm="AuxICA1")
out.AuxICA2 <- ICA(X, J=J, algorithm="AuxICA2")
out.IPCA <- ICA(X, J=J, algorithm="IPCA")
out.SIMBEC <- ICA(X, J=J, algorithm="SIMBEC")
out.AMUSE <- ICA(X, J=J, algorithm="AMUSE")
out.SOBI <- ICA(X, J=J, algorithm="SOBI")
out.FOBI <- ICA(X, J=J, algorithm="FOBI")
out.ProDenICA <- ICA(X, J=J, algorithm="ProDenICA")
out.RICA <- ICA(X, J=J, algorithm="RICA")

## Test Input object / type
### Test I-1: Object Names

expect_identical(names(formals(ICA)), c("X", "J", "algorithm", "num.iter", "thr", "verbose"))

### Test I-2: X

expect_identical(as.character(formals(ICA)$X), "")

### Test I-3: J

expect_identical(as.character(formals(ICA)$J), "")

### Test I-4: algorithm

expect_identical(formals(ICA)$algorithm,
    c("JADE", "AuxICA1", "AuxICA2", "IPCA", "SIMBEC", "AMUSE",
        "SOBI", "FOBI", "ProDenICA", "RICA"))

### Test I-5: num.iter

expect_identical(formals(ICA)$num.iter, 30)

### Test I-6: thr

expect_identical(formals(ICA)$thr, 1E-10)

### Test I-7: verbose

expect_identical(formals(ICA)$verbose, FALSE)

## Test Output object / type
### Test O-1: Object

expect_identical(is.list(out.JADE), TRUE)
expect_identical(is.list(out.AuxICA1), TRUE)
expect_identical(is.list(out.AuxICA2), TRUE)
expect_identical(is.list(out.IPCA), TRUE)
expect_identical(is.list(out.SIMBEC), TRUE)
expect_identical(is.list(out.AMUSE), TRUE)
expect_identical(is.list(out.SOBI), TRUE)
expect_identical(is.list(out.FOBI), TRUE)
expect_identical(is.list(out.ProDenICA), TRUE)
expect_identical(is.list(out.RICA), TRUE)

### Test O-2: Object Names

expect_identical(names(out.JADE),
    c("A", "S", "J", "algorithm", "num.iter", "thr", "verbose", "RecError", "RelChange"))
expect_identical(names(out.AuxICA1),
    c("A", "S", "J", "algorithm", "num.iter", "thr", "verbose", "RecError", "RelChange"))
expect_identical(names(out.AuxICA2),
    c("A", "S", "J", "algorithm", "num.iter", "thr", "verbose", "RecError", "RelChange"))
expect_identical(names(out.IPCA),
    c("A", "S", "J", "algorithm", "num.iter", "thr", "verbose", "RecError", "RelChange"))
expect_identical(names(out.SIMBEC),
    c("A", "S", "J", "algorithm", "num.iter", "thr", "verbose", "RecError", "RelChange"))
expect_identical(names(out.AMUSE),
    c("A", "S", "J", "algorithm", "num.iter", "thr", "verbose", "RecError", "RelChange"))
expect_identical(names(out.SOBI),
    c("A", "S", "J", "algorithm", "num.iter", "thr", "verbose", "RecError", "RelChange"))
expect_identical(names(out.FOBI),
    c("A", "S", "J", "algorithm", "num.iter", "thr", "verbose", "RecError", "RelChange"))
expect_identical(names(out.ProDenICA),
    c("A", "S", "J", "algorithm", "num.iter", "thr", "verbose", "RecError", "RelChange"))
expect_identical(names(out.RICA),
    c("A", "S", "J", "algorithm", "num.iter", "thr", "verbose", "RecError", "RelChange"))

### Test O-3: A

expect_identical(is.matrix(out.JADE$A), TRUE)
expect_identical(is.matrix(out.AuxICA1$A), TRUE)
expect_identical(is.matrix(out.AuxICA2$A), TRUE)
expect_identical(is.matrix(out.IPCA$A), TRUE)
expect_identical(is.matrix(out.SIMBEC$A), TRUE)
expect_identical(is.matrix(out.AMUSE$A), TRUE)
expect_identical(is.matrix(out.SOBI$A), TRUE)
expect_identical(is.matrix(out.FOBI$A), TRUE)
expect_identical(is.matrix(out.ProDenICA$A), TRUE)
expect_identical(is.matrix(out.RICA$A), TRUE)

expect_identical(dim(out.JADE$A), c(nrow(X), J))
expect_identical(dim(out.AuxICA1$A), c(nrow(X), J))
expect_identical(dim(out.AuxICA2$A), c(nrow(X), J))
expect_identical(dim(out.IPCA$A), c(nrow(X), J))
expect_identical(dim(out.SIMBEC$A), c(nrow(X), J))
expect_identical(dim(out.AMUSE$A), c(nrow(X), J))
expect_identical(dim(out.SOBI$A), c(nrow(X), J))
expect_identical(dim(out.FOBI$A), c(nrow(X), J))
expect_identical(dim(out.ProDenICA$A), c(nrow(X), J))
expect_identical(dim(out.RICA$A), c(nrow(X), J))

### Test O-4: S

expect_identical(is.matrix(out.JADE$S), TRUE)
expect_identical(is.matrix(out.AuxICA1$S), TRUE)
expect_identical(is.matrix(out.AuxICA2$S), TRUE)
expect_identical(is.matrix(out.IPCA$S), TRUE)
expect_identical(is.matrix(out.SIMBEC$S), TRUE)
expect_identical(is.matrix(out.AMUSE$S), TRUE)
expect_identical(is.matrix(out.SOBI$S), TRUE)
expect_identical(is.matrix(out.FOBI$S), TRUE)
expect_identical(is.matrix(out.ProDenICA$S), TRUE)
expect_identical(is.matrix(out.RICA$S), TRUE)

expect_identical(dim(out.JADE$S), c(J, ncol(X)))
expect_identical(dim(out.AuxICA1$S), c(J, ncol(X)))
expect_identical(dim(out.AuxICA2$S), c(J, ncol(X)))
expect_identical(dim(out.IPCA$S), c(J, ncol(X)))
expect_identical(dim(out.SIMBEC$S), c(J, ncol(X)))
expect_identical(dim(out.AMUSE$S), c(J, ncol(X)))
expect_identical(dim(out.SOBI$S), c(J, ncol(X)))
expect_identical(dim(out.FOBI$S), c(J, ncol(X)))
expect_identical(dim(out.ProDenICA$S), c(J, ncol(X)))
expect_identical(dim(out.RICA$S), c(J, ncol(X)))


### Test O-5: J

expect_identical(out.JADE$J, J)
expect_identical(out.AuxICA1$J, J)
expect_identical(out.AuxICA2$J, J)
expect_identical(out.IPCA$J, J)
expect_identical(out.SIMBEC$J, J)
expect_identical(out.AMUSE$J, J)
expect_identical(out.SOBI$J, J)
expect_identical(out.FOBI$J, J)
expect_identical(out.ProDenICA$J, J)
expect_identical(out.RICA$J, J)

### Test O-6: algorithm

expect_identical(out.JADE$algorithm, "JADE")
expect_identical(out.AuxICA1$algorithm, "AuxICA1")
expect_identical(out.AuxICA2$algorithm, "AuxICA2")
expect_identical(out.IPCA$algorithm, "IPCA")
expect_identical(out.SIMBEC$algorithm, "SIMBEC")
expect_identical(out.AMUSE$algorithm, "AMUSE")
expect_identical(out.SOBI$algorithm, "SOBI")
expect_identical(out.FOBI$algorithm, "FOBI")
expect_identical(out.ProDenICA$algorithm, "ProDenICA")
expect_identical(out.RICA$algorithm, "RICA")

### Test O-7: num.iter

expect_identical(out.JADE$num.iter, formals(ICA)$num.iter)
expect_identical(out.AuxICA1$num.iter, formals(ICA)$num.iter)
expect_identical(out.AuxICA2$num.iter, formals(ICA)$num.iter)
expect_identical(out.IPCA$num.iter, formals(ICA)$num.iter)
expect_identical(out.SIMBEC$num.iter, formals(ICA)$num.iter)
expect_identical(out.AMUSE$num.iter, formals(ICA)$num.iter)
expect_identical(out.SOBI$num.iter, formals(ICA)$num.iter)
expect_identical(out.FOBI$num.iter, formals(ICA)$num.iter)
expect_identical(out.ProDenICA$num.iter, formals(ICA)$num.iter)
expect_identical(out.RICA$num.iter, formals(ICA)$num.iter)

### Test O-8: thr

expect_identical(out.JADE$thr, formals(ICA)$thr)
expect_identical(out.AuxICA1$thr, formals(ICA)$thr)
expect_identical(out.AuxICA2$thr, formals(ICA)$thr)
expect_identical(out.IPCA$thr, formals(ICA)$thr)
expect_identical(out.SIMBEC$thr, formals(ICA)$thr)
expect_identical(out.AMUSE$thr, formals(ICA)$thr)
expect_identical(out.SOBI$thr, formals(ICA)$thr)
expect_identical(out.FOBI$thr, formals(ICA)$thr)
expect_identical(out.ProDenICA$thr, formals(ICA)$thr)
expect_identical(out.RICA$thr, formals(ICA)$thr)


### Test O-9: verbose

expect_identical(out.JADE$verbose, formals(ICA)$verbose)
expect_identical(out.AuxICA1$verbose, formals(ICA)$verbose)
expect_identical(out.AuxICA2$verbose, formals(ICA)$verbose)
expect_identical(out.IPCA$verbose, formals(ICA)$verbose)
expect_identical(out.SIMBEC$verbose, formals(ICA)$verbose)
expect_identical(out.AMUSE$verbose, formals(ICA)$verbose)
expect_identical(out.SOBI$verbose, formals(ICA)$verbose)
expect_identical(out.FOBI$verbose, formals(ICA)$verbose)
expect_identical(out.ProDenICA$verbose, formals(ICA)$verbose)
expect_identical(out.RICA$verbose, formals(ICA)$verbose)

### Test O-10: RecError

expect_identical(is.vector(out.JADE$RecError), TRUE)
expect_identical(is.vector(out.AuxICA1$RecError), TRUE)
expect_identical(is.vector(out.AuxICA2$RecError), TRUE)
expect_identical(is.vector(out.IPCA$RecError), TRUE)
expect_identical(is.vector(out.SIMBEC$RecError), TRUE)
expect_identical(is.vector(out.AMUSE$RecError), TRUE)
expect_identical(is.vector(out.SOBI$RecError), TRUE)
expect_identical(is.vector(out.FOBI$RecError), TRUE)
expect_identical(is.vector(out.ProDenICA$RecError), TRUE)
expect_identical(is.vector(out.RICA$RecError), TRUE)


### Test O-11: RelChange

expect_identical(is.vector(out.JADE$RelChange), TRUE)
expect_identical(is.vector(out.AuxICA1$RelChange), TRUE)
expect_identical(is.vector(out.AuxICA2$RelChange), TRUE)
expect_identical(is.vector(out.IPCA$RelChange), TRUE)
expect_identical(is.vector(out.SIMBEC$RelChange), TRUE)
expect_identical(is.vector(out.AMUSE$RelChange), TRUE)
expect_identical(is.vector(out.SOBI$RelChange), TRUE)
expect_identical(is.vector(out.FOBI$RelChange), TRUE)
expect_identical(is.vector(out.ProDenICA$RelChange), TRUE)
expect_identical(is.vector(out.RICA$RelChange), TRUE)


## Test Error
### Test E-1: X

expect_error(ICA(as.data.frame(X), J=J))

### Test E-2: J

expect_error(ICA(X, J="5"))
expect_error(ICA(X, J=c(2,4)))
expect_error(ICA(X, J=10^10)

### Test E-3: num.iter

expect_error(ICA(X, J=J, num.iter="100"))
expect_error(ICA(X, J=J, num.iter=-1))

### Test E-4: thr

expect_error(ICA(X, J=J, thr="0.1"))

### Test E-5: verbose

expect_error(ICA(X, J=J, verbose="verbose"))

## Test Decrease of Error

### Test D-1: RecError
.sampleRank <- function(x){
	rank(c(x[2], median(x), rev(x)[1]))
}

expect_identical(.sampleRank(out.JADE$RecError), 3:1)
expect_identical(.sampleRank(out.AuxICA1$RecError), 3:1)
expect_identical(.sampleRank(out.AuxICA2$RecError), 3:1)
expect_identical(.sampleRank(out.IPCA$RecError), 3:1)
expect_identical(.sampleRank(out.SIMBEC$RecError), 3:1)
expect_identical(.sampleRank(out.AMUSE$RecError), 3:1)
expect_identical(.sampleRank(out.SOBI$RecError), 3:1)
expect_identical(.sampleRank(out.FOBI$RecError), 3:1)
expect_identical(.sampleRank(out.ProDenICA$RecError), 3:1)
expect_identical(.sampleRank(out.RICA$RecError), 3:1)

### Test D-2: RelChange

expect_identical(.sampleRank(out.JADE$RelChange), 3:1)
expect_identical(.sampleRank(out.AuxICA1$RelChange), 3:1)
expect_identical(.sampleRank(out.AuxICA2$RelChange), 3:1)
expect_identical(.sampleRank(out.IPCA$RelChange), 3:1)
expect_identical(.sampleRank(out.SIMBEC$RelChange), 3:1)
expect_identical(.sampleRank(out.AMUSE$RelChange), 3:1)
expect_identical(.sampleRank(out.SOBI$RelChange), 3:1)
expect_identical(.sampleRank(out.FOBI$RelChange), 3:1)
expect_identical(.sampleRank(out.ProDenICA$RelChange), 3:1)
expect_identical(.sampleRank(out.RICA$RelChange), 3:1)

## Session Information

sessionInfo()
