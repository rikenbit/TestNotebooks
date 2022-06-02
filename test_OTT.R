# Test code of TOT Function

## Dependent Packages Installation

if (!requireNamespace("testthat", quietly = TRUE)){
    install.packages("testthat", repos="https://cran.r-project.org")
}
library("testthat")
if (!requireNamespace("rTensor", quietly = TRUE)){
    install.packages("rTensor", repos="https://cran.r-project.org")
}
library("rTensor")

## OTT Function (Paste your OTT code here)

# num.sample / num.iterは論文の値にしましたが、大きすぎたら、小さめの値に変えてください
OTT <- function(X, Y, ps=NULL, qs=NULL, loss=dist, num.sample=1000,
    num.iter=1000, thr=1e-10, epsilon=1e-10, verbose=FALSE){
    ######################################
    # Argument Check
    ######################################
    .checkOTT(X, Y, ps, qs, loss, num.sample, num.iter, thr, epsilon, verbose)
    ######################################
    # Initialization
    ######################################
    int <- .initOTT(X, Y, ps, qs, thr)
    Ls <- int$Ls
    ps <- int$ps
    qs <- int$qs
    Ts <- int$Ts
    A <- int$A
    RecError <- int$RecError
    RelChange <- int$RelChange
    ######################################
    # Iteration
    ######################################
    iter <- 1
    while ((RelChange[iter] > thr) && (iter <= num.iter)) {
        for(a in seq(A)){
            Ts[[a]] <- ...
         }
        iter <- iter + 1
        # 目的関数の値を入れる
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
    list(Ts=Ts, num.sample=num.sample, num.iter=num.iter,
        thr=thr, epsilon=epsilon, verbose=verbose,
        RecError=RecError, RelChange=RelChange)
}

.checkOTT <- function(X, Y, ps, qs, loss, num.sample,
    num.iter, thr, epsilon, verbose){
    stopifnot(is(X)[1] == "Tensor")
    stopifnot(is(Y)[1] == "Tensor")
    stopifnot(identical(dim(X), dim(Y)))
    if(!is.null(ps)){
        stopifnot(is.list(ps))
        l1 <- dim(X)
        l2 <- lapply(ps, length)
        stopifnot(identical(l1, l2))
    }
    if(!is.null(qs)){
        stopifnot(is.list(qs))
        l1 <- dim(Y)
        l2 <- lapply(qs, length)
        stopifnot(identical(l1, l2))
    }
    stopifnot(is(loss)[1] == "function")
    stopifnot(is.numeric(num.sample))
    stopifnot(num.sample >= 1)
    stopifnot(num.sample <= ...) # 上限値をここに入れてください
    stopifnot(is.numeric(num.iter))
    stopifnot(num.iter >= 1)
    stopifnot(is.numeric(thr))
    stopifnot(thr >= 0)
    stopifnot(is.numeric(epsilon))
    stopifnot(epsilon >= 0)
    stopifnot(is.logical(verbose))
}

.initOTT <- function(X, Y, ps, qs, thr){
    # Lossの計算
    Ls <- ...
    # 最適輸送計画
    if (is.null(ps)){
        # 何も指定がない場合は、一様分布
        ps <- lapply(dim(X), function(x){
            out <- rep(1, length=x)
            out / sum(out)
        })
    }
    if (is.null(qs)){
        # 何も指定がない場合は、一様分布
        qs <- lapply(dim(Y), function(x){
            out <- rep(1, length=x)
            out / sum(out)
        })
    }
    A <- length(dim(X))
    Ts <- lapply(seq(A), function(a){
        ps[[a]] %o% qs[[a]]
    })
    # Reconstruction Error
    RecError <- ...
    # Relative Change
    RelChange <- thr * 10
    list(Ls=Ls, ps=ps, qs=qs, Ts=Ts, A=A,
        RecError=RecError, RelChange=RelChange)
}


## Simulation Datasets

arrX <- array(runif(10*12*14), dim=c(10,12,14))
arrY <- array(array(runif(11*13*15), dim=c(11,13,15)))
X <- as.tensor(arrX)
Y <- as.tensor(arrY)

## Perform OTT against Simulation Dataset

out <- OTT(X, Y)


## Test Input object / type
### Test I-1: Object Names

expect_identical(names(formals(OTT)),
    c("X", "Y", "ps", "qs", "loss", "num.sample",
        "num.iter", "thr", "epsilon", "verbose"))

### Test I-2: X

expect_identical(as.character(formals(OTT)$X), "")

### Test I-3: Y

expect_identical(as.character(formals(OTT)$Y), "")

### Test I-4: ps

expect_identical(formals(OTT)$ps, NULL)

### Test I-5: qs

expect_identical(formals(OTT)$qs, NULL)

### Test I-5: loss

expect_true(is.function(formals(OTT)$loss))

### Test I-6: num.sample

expect_identical(formals(OTT)$num.sample, 1000)

### Test I-7: num.iter

expect_identical(formals(OTT)$num.iter, 1000)

### Test I-8: thr

expect_identical(formals(OTT)$thr, 1E-10)

### Test I-9: epsilon

expect_identical(formals(OTT)$epsilon, 1E-10)

### Test I-10: verbose

expect_identical(formals(OTT)$verbose, FALSE)


## Test Output object / type
### Test O-1: Object

expect_identical(is.list(out), TRUE)

### Test O-2: Object Names

expect_identical(names(out),
    c("Ts", "num.sample", "num.iter", "thr",
        "epsilon", "verbose", "RecError", "RelChange"))

### Test 0-3: Ts

expect_identical(is.list(out$Ts), TRUE)
expect_identical(is.matrix(out$Ts[[1]]), TRUE)
expect_identical(is.matrix(out$Ts[[2]]), TRUE)
expect_identical(is.matrix(out$Ts[[3]]), TRUE)

expect_identical(dim(out$Ts[[1]]), c(dim(X)[1], dim(Y)[1]))
expect_identical(dim(out$Ts[[2]]), c(dim(X)[2], dim(Y)[2]))
expect_identical(dim(out$Ts[[3]]), c(dim(X)[3], dim(Y)[3]))

### Test 0-4: num.sample

expect_identical(out$num.sample, formals(OTT)$num.sample)

### Test 0-5: num.iter

expect_identical(out$num.iter, formals(OTT)$num.iter)

### Test 0-6: thr

expect_identical(out$thr, formals(OTT)$thr)

### Test 0-7: epsilon

expect_identical(out$epsilon, formals(OTT)$epsilon)

### Test 0-8: verbose

expect_identical(out$verbose, formals(OTT)$verbose)

### Test 0-9: RecError

expect_identical(is.vector(out$RecError), TRUE)

### Test 0-10: RelChange

expect_identical(is.vector(out$RelChange), TRUE)


## Test Error
### Test E-1: X / Y

expect_error(OTT(arrX, Y))
expect_error(OTT(X, arrY))
expect_error(OTT(arrX, arrY))

### Test E-2: ps

expect_error(OTT(X, Y, ps=runif(1:10)))

### Test E-3: qs

expect_error(OTT(X, Y, qs=runif(1:10)))


### Test E-4: loss

expect_error(OTT(X, Y, loss="loss"))
expect_error(OTT(X, Y, loss=marix(runif(3*4), nrow=3, ncol=4)))

### Test E-5: num.sample

expect_error(OTT(X, Y, num.sample="100"))
expect_error(OTT(X, Y, num.sample=-1))

### Test E-6: num.iter

expect_error(OTT(X, Y, num.iter="100"))
expect_error(OTT(X, Y, num.iter=-1))

### Test E-7: thr

expect_error(OTT(X, Y, thr="100"))
expect_error(OTT(X, Y, epsilon=-1))

### Test E-8: epsilon

expect_error(OTT(X, Y, epsilon="100"))
expect_error(OTT(X, Y, epsilon=-1))

### Test E-9: verbose

expect_error(OTT(X, Y, verbose="verbose"))




## Test Decrease of Error
### Test D-1: RecError
.sampleRank <- function(x){
	rank(c(x[2], median(x), rev(x)[1]))
}
expect_identical(.sampleRank(out$RecError), 3:1)



### Test D-2: RelChange
expect_identical(.sampleRank(out$RelChange), 3:1)


## Session Information
sessionInfo()
