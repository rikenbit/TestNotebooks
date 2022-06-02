# Test code of GroupICA Function

## Dependent Packages Installation
if (!requireNamespace("testthat", quietly = TRUE)){
    install.packages("testthat", repos="https://cran.r-project.org")
}
library("testthat")
if (!requireNamespace("rlist", quietly = TRUE)){
    install.packages("rlist", repos="https://cran.r-project.org")
}
library("rlist")

## GroupICA Function (Paste your GroupICA code here)

GroupICA <- function(Xs, J1, J2=J1,
    algorithm=c("pooled", "Calhourn2009", "Pfister2018"),
    ica.algorithm=c("JADE", "AuxICA1", "AuxICA2", "IPCA",
        "SIMBEC", "AMUSE",
        "SOBI", "FOBI", "ProDenICA", "RICA"),
    num.iter=30, thr=1E-10, verbose=FALSE){
    ######################################
    # Argument Check
    ######################################
    .checkGroupICA(Xs, J1, J2, algorithm, num.iter, thr, verbose)
    algorithm <- match.arg(algorithm)
    ica.algorithm <- match.arg(algorithm)
    l <- length(Xs)
    if(algorithm == "pooled"){
        Xpooled <- rlist::list.cbind(Xs)
        out <- ICA(Xpooled, J=J1, algorithm=ica.algorithm,
            num.iter=num.iter, thr=thr, verbose=verbose)
        # Output
        A <- out$A
        Ss <- lapply(seq(l), function(x){out$S})
        RecError <- out$RecError
        RelChange <- out$RelChange
    }
    if(algorithm == "Calhourn2009"){
        # 1. Each PCA
        res.pca <- lapply(Xs, function(x){
            prcomp(x, center=TRUE, scale=FALSE)
        })
        PCcat <- rlist::list.cbind(lapply(res.pca, function(x){x$x[, seq(J1)]}))
        # 2. Merged PCA
        res.merged.pca <- prcomp(PCcat, center=TRUE, scale=FALSE)
        # 3. Merged ICA
        res.merged.ica <- ICA(res.merged.pca$x[, seq(J2)],
            J=J2, algorithm=ica.algorithm,
            num.iter=num.iter, thr=thr, verbose=verbose)
        # 4. Reconstruction / Output
        A <- res.merged.ica$A
        Ss <- lapply(seq(l), function(x){
            idx <- .eachidx(J1, l, x)
            res.merged.ica$S %*%
                t(res.merged.pca$rotation[idx, seq(J2)]) %*%
                    t(res.pca[[x]]$rotation[, seq(J1)])
        })
        RecError <- NULL
        RelChange <- NULL
    }
    if(algorithm == "Pfister2018"){      
        ######################################
        # Initialization (e.g. Whiteniing)
        ######################################
        int <- .initGroupICA(Xs)
        X <- int$X # poolされて出力されるように
        g <- int$g # group index
        Pg <- int$Pg # group-wise partition（GroupICAのオプションにする必要有り）
        M <- int$M # empty list
        for(g in G){
            for(e in E){
                M[...] <- cov(Xe) - cov(Xe)
            }
        }
        out <- .ApproximateJointDiagonalizer(M)
        # Output
        A <- X %*% t(out)
        Ss <- out # 個体ごとに要素を持つリストにする
        RecError <- NULL
        RelChange <- NULL
    }    
    # Output
    list(A=A, Ss=Ss, J1=J1, J2=J2,
        algorithm=algorithm, ica.algorithm=ica.algorithm,
        num.iter=num.iter, thr=thr, verbose=verbose,
        RecError=RecError, RelChange=RelChange)
}

.checkGroupICA <- function(Xs, J1, J2, algorithm, num.iter, thr, verbose){
    stopifnot(is.list(Xs))
    nr <- lapply(Xs, nrow)
    all.equal(length(unique(nr)), 1)
    stopifnot(is.numeric(J1))
    stopifnot(length(J1) == 1)
    lapply(Xs, function(x){
        stopifnot(min(dim(x)) >= J1)
    })
    if(algorithm == "Calhourn2009"){
        stopifnot(is.numeric(J2))
        stopifnot(length(J2) == 1)
        stopifnot(min(nrow(Xs[[1]]), length(Xs)*J1) >= J2)
    }
    stopifnot(is.numeric(num.iter))
    stopifnot(num.iter > 0)
    stopifnot(is.numeric(thr))
    stopifnot(is.logical(verbose))
}

.initGroupICA <- function(Xs){
    X <- ...
    g <- ...
    Pg <- ...
    M <- ...
    list(X=X, g=g, Pg=Pg, M=M)
}

.eachidx <- function(J1, l, x){
    out <- 1:(J1*l)
    start <- seq(from=1, to=J1*l, by=J1)[x]
    end <- seq(from=J1, to=J1*l, by=J1)[x]
    out[start:end]
}

.ApproximateJointDiagonalizer <- function(M){
    ...
}



## Simulation Dataset

X1 <- matrix(runif(100*200), nrow=100, ncol=200)
X2 <- matrix(runif(100*150), nrow=100, ncol=150)
X3 <- matrix(runif(100*250), nrow=100, ncol=250)

Xs <- list(X1=X1, X2=X2, X3=X3)


## Perform GroupICA against Simulation Dataset

J1 <- 5
out.pooled <- GroupICA(Xs, J1=J1, algorithm="pooled")
out.Calhourn2009 <- GroupICA(Xs, J1=J1, algorithm="Calhourn2009")
out.Pfister2018 <- GroupICA(Xs, J1=J1, algorithm="Pfister2018")

## Test Input object / type
### Test I-1: Object Names

expect_identical(names(formals(GroupICA)),
    c("Xs", "J1", "J2", "algorithm", "ica.algorithm", "num.iter", "thr", "verbose"))

### Test I-2: Xs

expect_identical(as.character(formals(GroupICA)$Xs), "")

### Test I-3: J1

expect_identical(as.character(formals(GroupICA)$J1), "")

### Test I-4: J2

expect_identical(as.character(formals(GroupICA)$J2), "")

### Test I-5: algorithm

expect_identical(formals(GroupICA)$algorithm,
    c("pooled", "Calhourn2009", "Pfister2018"))

### Test I-6: ica.algorithm

expect_identical(formals(GroupICA)$ica.algorithm,
    c("JADE", "AuxICA1", "AuxICA2", "IPCA",]
        "SIMBEC", "AMUSE", "SOBI", "FOBI", "ProDenICA", "RICA"))

### Test I-7: num.iter

expect_identical(formals(GroupICA)$num.iter, 30)

### Test I-8: thr

expect_identical(formals(GroupICA)$thr, 1E-10)

### Test I-9: verbose

expect_identical(formals(GroupICA)$verbose, FALSE)


## Test Output object / type
### Test O-1: Object

expect_identical(is.list(out.pooled), TRUE)
expect_identical(is.list(out.Calhourn2009), TRUE)
expect_identical(is.list(out.Pfister2018), TRUE)

### Test O-2: Object Names

expect_identical(names(out.pooled),
    c("A", "Ss", "J1", "J2", "algorithm", "ica.algorithm",
        "num.iter", "thr", "verbose", "RecError", "RelChange"))
expect_identical(names(out.Calhourn2009),
    c("A", "Ss", "J1", "J2", "algorithm", "ica.algorithm",
        "num.iter", "thr", "verbose", "RecError", "RelChange"))
expect_identical(names(out.Pfister2018),
    c("A", "Ss", "J1", "J2", "algorithm", "ica.algorithm",
        "num.iter", "thr", "verbose", "RecError", "RelChange"))

### Test O-3: A

expect_identical(is.matrix(out.pooled$A), TRUE)
expect_identical(is.matrix(out.Calhourn2009$A), TRUE)
expect_identical(is.matrix(out.Pfister2018$A), TRUE)

expect_identical(dim(out.pooled$A), c(nrow(X[[1]]), J1))
expect_identical(dim(out.Calhourn2009$A), c(nrow(X[[1]]), J1))
expect_identical(dim(out.Pfister2018$A), c(nrow(X[[1]]), J1))

### Test O-4: Ss

expect_identical(is.matrix(out.pooled$Ss[[1]]), TRUE)
expect_identical(is.matrix(out.pooled$Ss[[2]]), TRUE)
expect_identical(is.matrix(out.pooled$Ss[[3]]), TRUE)
expect_identical(is.matrix(out.Calhourn2009$Ss[[1]]), TRUE)
expect_identical(is.matrix(out.Calhourn2009$Ss[[2]]), TRUE)
expect_identical(is.matrix(out.Calhourn2009$Ss[[3]]), TRUE)
expect_identical(is.matrix(out.Pfister2018$Ss[[1]]), TRUE)
expect_identical(is.matrix(out.Pfister2018$Ss[[2]]), TRUE)
expect_identical(is.matrix(out.Pfister2018$Ss[[3]]), TRUE)

expect_identical(dim(out.pooled$Ss[[1]]), c(J1, ncol(X[[1]])))
expect_identical(dim(out.pooled$Ss[[2]]), c(J1, ncol(X[[2]])))
expect_identical(dim(out.pooled$Ss[[3]]), c(J1, ncol(X[[3]])))
expect_identical(dim(out.Calhourn2009$Ss[[1]]), c(J1, ncol(X[[1]])))
expect_identical(dim(out.Calhourn2009$Ss[[2]]), c(J1, ncol(X[[2]])))
expect_identical(dim(out.Calhourn2009$Ss[[3]]), c(J1, ncol(X[[3]])))
expect_identical(dim(out.Pfister2018$Ss[[1]]), c(J1, ncol(X[[1]])))
expect_identical(dim(out.Pfister2018$Ss[[2]]), c(J1, ncol(X[[2]])))
expect_identical(dim(out.Pfister2018$Ss[[3]]), c(J1, ncol(X[[3]])))

### Test O-5: J1

expect_identical(out.pooled$J1, J1)
expect_identical(out.Calhourn2009$J1, J1)
expect_identical(out.Pfister2018$J1, J1)

### Test O-6: J2

expect_identical(out.pooled$J2, J1)
expect_identical(out.Calhourn2009$J2, J1)
expect_identical(out.Pfister2018$J2, J1)

### Test O-7: algorithm

expect_identical(out.pooled$algorithm, "pooled")
expect_identical(out.Calhourn2009$algorithm, "Calhourn2009")
expect_identical(out.Pfister2018$algorithm, "Pfister2018")

### Test O-8: ica.algorithm

expect_identical(out.pooled$ica.algorithm, "JADE")
expect_identical(out.Calhourn2009$ica.algorithm, "JADE")
expect_identical(out.Pfister2018$ica.algorithm, "JADE")

### Test O-9: num.iter

expect_identical(out.pooled$num.iter, formals(GroupICA)$num.iter)
expect_identical(out.Calhourn2009$num.iter, formals(GroupICA)$num.iter)
expect_identical(out.Pfister2018$num.iter, formals(GroupICA)$num.iter)

### Test O-10: thr

expect_identical(out.pooled$thr, formals(GroupICA)$thr)
expect_identical(out.Calhourn2009$thr, formals(GroupICA)$thr)
expect_identical(out.Pfister2018$thr, formals(GroupICA)$thr)

### Test O-11: verbose

expect_identical(out.pooled$verbose, formals(GroupICA)$verbose)
expect_identical(out.Calhourn2009$verbose, formals(GroupICA)$verbose)
expect_identical(out.Pfister2018$verbose, formals(GroupICA)$verbose)

### Test O-12: RecError

expect_identical(is.vector(out.pooled$RecError), TRUE)
expect_identical(is.null(out.Calhourn2009$RecError), TRUE)
expect_identical(is.null(out.Pfister2018$RecError), TRUE)

### Test O-13: RelChange

expect_identical(is.vector(out.pooled$RelChange), TRUE)
expect_identical(is.null(out.Calhourn2009$RelChange), TRUE)
expect_identical(is.null(out.Pfister2018$RelChange), TRUE)

## Test Error
### Test E-1: X

expect_error(GroupICA(X, J1=J1))

### Test E-2: J1

expect_error(GroupICA(Xs, J1="5"))
expect_error(GroupICA(Xs, J1=c(2,4)))
expect_error(GroupICA(Xs, J1=10^10)

### Test E-3: J2

expect_error(GroupICA(Xs, J2="5"))
expect_error(GroupICA(Xs, J2=c(2,4)))
expect_error(GroupICA(Xs, J2=10^10)

### Test E-4: algorithm

expect_error(GroupICA(Xs, algorithm="poooled"))

### Test E-5: ica.algorithm

expect_error(GroupICA(Xs, ica.algorithm="JAAE"))

### Test E-6: num.iter

expect_error(GroupICA(Xs, J=J, num.iter="100"))
expect_error(GroupICA(Xs, J=J, num.iter=-1))

### Test E-7: thr

expect_error(GroupICA(X, J=J, thr="0.1"))

### Test E-8: verbose

expect_error(GroupICA(X, J=J, verbose="verbose"))


## Test Decrease of Error

### Test D-1: RecError
.sampleRank <- function(x){
    rank(c(x[2], median(x), rev(x)[1]))
}

expect_identical(.sampleRank(out.pooled$RecError), 3:1)

### Test D-2: RelChange

expect_identical(.sampleRank(out.pooled$RelChange), 3:1)


## Session Information

sessionInfo()
