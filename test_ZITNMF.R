# Test code of ZITNMF Function

## Dependent Packages Installation
if (!requireNamespace("testthat", quietly = TRUE)){
    install.packages("testthat", repos="https://cran.r-project.org")
}
library("testthat")
if (!requireNamespace("nnTensor", quietly = TRUE)){
    install.packages("nnTensor", repos="https://cran.r-project.org")
}
library("nnTensor")

## ZITNMF Function (Paste your ZITNMF code here)

ZITNMF <- function(X, Z=NULL, pseudocount=1e-10,
	initF=NULL, initA=NULL, fixF=FALSE, fixA=FALSE,
    init = c("NMF", "ALS", "Random"),
	J=3, Beta=2, phi=1, thr=1e-10, num.iter=100, verbose=FALSE){
    ######################################
    # Argument Check
    ######################################
    .checkZITNMF(X, Z, pseudocount, initF, initA, fixF, fixA,
        J, Beta, phi, thr, num.iter, verbose)
    init <- match.arg(init)
    ######################################
    # Initialization
    ######################################
    int <- .initZITNMF(X, Z, pseudocount, initF, initA, J, Beta)
    X <- int$X
    E <- int$E
    F <- int$F
    A <- int$A
    w <- int$w
    RecError <- int$RecError
    RelChange <- int$RelChange
    ######################################
    # Iteration
    ######################################
	iter <- 1
    while ((RecError[iter] > thr) && (iter <= num.iter)) {
    	# Update Z
    	Z <- ...
    	# Update w
    	w <- mean(Z)
    	if(!fixF){
			numerF <- ...
			denomF <- ...
			F <- F * (numerF / denomF)^.rho(Beta)
    	}
    	if(!fixA){
			numerA <- ...
			denomA <- ...
			A <- A * (numerA / denomA)^.rho(Beta)
    	}
        # After Update U, V
        iter <- iter + 1
        # ここは、LogLikelihoodの方が良いかもしれません
        # （.recErrorは最小自乗誤差しか計算していないので）
        X_bar <- nnTensor:::.recMatrix(F, A)
        RecError[iter] <- nnTensor:::.recError(X, X_bar)
        TrainRecError[iter] <- nnTensor:::.recError(Z*X, Z*X_bar)
        TestRecError[iter] <- nnTensor:::.recError((1-Z)*X, (1-Z)*X_bar)
        RelChange[iter] <- abs(pre_Error - RecError[iter]) / RecError[iter]
        if (verbose) {
            cat(paste0(iter-1, " / ", num.iter, " |Previous Error - Error| / Error = ",
                RelChange[iter], "\n"))
        }
        if (is.nan(RelChange[iter])) {
            stop("NaN is generated. Please run again or change the parameters.\n")
        }
    }
    names(RecError) <- c("offset", 1:(iter-1))
    names(TrainRecError) <- c("offset", 1:(iter-1))
    names(TestRecError) <- c("offset", 1:(iter-1))
    names(RelChange) <- c("offset", 1:(iter-1))
    list(F=F, A=A, Z=Z, w=w, RecError=RecError, TrainRecError=TrainRecError,
    	TestRecError=TestsRecError, RelChange=RelChange)
}


.checkZITNMF <- function(X, Z, pseudocount, initF, initA, fixF, fixA,
    J, Beta, phi, thr, num.iter, verbose){
    stopifnot(is.matrix(X))
    if(!is.null(Z)){
    	stopifnot(identical(dim(X), dim(Z)))
    }
    stopifnot(is.numeric(pseudocount))
    stopifnot(pseudocount >= 0)
    if(!is.null(initF)){
    	stopifnot(is.matrix(initF))
    	stopifnot(dim(initF)[1] == nrow(X))
    	stopifnot(dim(initF)[2] == J)
    }
    if(!is.null(initA)){
    	stopifnot(is.matrix(initA))
    	stopifnot(dim(initA)[1] == J)
    	stopifnot(dim(initA)[2] == ncol(X))
    }
    stopifnot(is.logical(fixF))
    stopifnot(is.logical(fixA))
    stopifnot(is.numeric(J))
    stopifnot(J <= min(dim(X)))
    stopifnot(is.numetic(Beta))
    stopifnot(is.numetic(phi))
    stopifnot(is.numeric(thr))
    stopifnot(thr > 0)
    stopifnot(is.numetic(num.iter))
    stopifnot(num.iter >= 0)
    stopifnot(is.logical(verbose))
}

.initZITNMF <- function(X, Z, pseudocount, initF, initA, J, Beta){
    X[which(X == 0)] <- pseudocount
	E <- X
	E[] <- 1
	if(is.null(initF) || is.null(initA)){
		if(init == "NMF"){
		    res.NMF <- NMF(X, J = J, algorithm = "Beta", Beta=Beta)
		    if(is.null(initF)){
				initF <- res.NMF$U
		    }
		    if(is.null(initA)){
				initA <- t(res.NMF$V)
		    }			
		}
		if(init == "ALS"){
			res.ALS <- svd(X)
		    if(is.null(initF)){
				initF <- nnTensor:::.positive(res.ALS$u[, seq(J)])
		    }
		    if(is.null(initA)){
				initA <- t(nnTensor:::.positive(res.ALS$v[, seq(J)]))
		    }						
		}
		if(init == "Random"){
		    if(is.null(initF)){
				initF <- matrix(runif(nrow(X)*J), nrow=nrow(X), ncol=J)
		    }
		    if(is.null(initA)){
				initA <- matrix(runif(J*ncol(X)), nrow=J, ncol=ncol(X))
		    }						
		}
	}
	w <- ...
	list(X=X, E=E, F=initF, A=initA, w=w)
}

.rho <- function(Beta){
    if(Beta < 1){
        rho_beta <- 1 / (2 - Beta)
    }
    if((1 <= Beta) && (Beta <= 2)){
        rho_beta <- 1
    }
    if(Beta > 2){
        rho_beta <- 1 / (Beta - 1)
    }
    rho_beta
}

## Simulation Dataset

X <- nnTensor::toyModel("NMF")

## Perform ZITNMF against Simulation Dataset

J <- 5
out.Frobenius <- ZITNMF(X, J=J, Beta=2)
out.KL <- ZITNMF(X, J=J, Beta=1)
out.IS <- ZITNMF(X, J=J, Beta=0)

## Test Input object / type
### Test I-1: Object Names

expect_identical(names(formals(ZITNMF)),
	c("X", "Z" ,"pseudocount", "initF", "initA", "fixF", "fixA",
		"init", "J", "Beta", "phi", "thr", "num.iter", "verbose"))

### Test I-2: X

expect_identical(as.character(formals(ZITNMF)$X), "")

### Test I-3: Z

expect_identical(formals(ZITNMF)$Z, NULL)

### Test I-4: pseudocount

expect_identical(formals(ZITNMF)$pseudocount, 1e-10)

### Test I-5: initF

expect_identical(formals(ZITNMF)$initF, NULL)

### Test I-6: initA

expect_identical(formals(ZITNMF)$initA, NULL)

### Test I-7: fixF

expect_identical(formals(ZITNMF)$fixF, FALSE)

### Test I-8: fixA

expect_identical(formals(ZITNMF)$fixA, FALSE)

### Test I-9: init

expect_identical(formals(ZITNMF)$init, c("NMF", "ALS", "Random"))

### Test I-10: J

expect_identical(formals(ZITNMF)$J, 3)

### Test I-11: Beta

expect_identical(formals(ZITNMF)$Beta, 2)

### Test I-12: phi

expect_identical(formals(ZITNMF)$phi, 1)

### Test I-13: thr

expect_identical(formals(ZITNMF)$thr, 1e-10)

### Test I-14: num.iter

expect_identical(formals(ZITNMF)$num.iter, 100)

### Test I-15: verbose

expect_identical(formals(ZITNMF)$verbose, FALSE)


## Test Output object / type
### Test O-1: Object

expect_identical(is.list(out.Frobenius), TRUE)
expect_identical(is.list(out.KL), TRUE)
expect_identical(is.list(out.IS), TRUE)

### Test O-2: Object Names

expect_identical(names(out.Frobenius),
    c("F", "A", "Z", "w", "RecError", "TrainRecError", "TestsRecError", "RelChange"))
expect_identical(names(out.KL),
    c("F", "A", "Z", "w", "RecError", "TrainRecError", "TestsRecError", "RelChange"))
expect_identical(names(out.IS),
    c("F", "A", "Z", "w", "RecError", "TrainRecError", "TestsRecError", "RelChange"))

### Test O-3: F

expect_identical(is.matrix(out.Frobenius$F), TRUE)
expect_identical(is.matrix(out.KL$F), TRUE)
expect_identical(is.matrix(out.IS$F), TRUE)

expect_identical(dim(out.Frobenius$F), c(nrow(X), J))
expect_identical(dim(out.KL$F), c(nrow(X), J))
expect_identical(dim(out.IS$F), c(nrow(X), J))

### Test O-4: A

expect_identical(is.matrix(out.Frobenius$A), TRUE)
expect_identical(is.matrix(out.KL$A), TRUE)
expect_identical(is.matrix(out.IS$A), TRUE)

expect_identical(dim(out.Frobenius$A), c(J, ncol(X)))
expect_identical(dim(out.KL$A), c(J, ncol(X)))
expect_identical(dim(out.IS$A), c(J, ncol(X)))

### Test O-5: Z

expect_identical(is.matrix(out.Frobenius$Z), TRUE)
expect_identical(is.matrix(out.KL$Z), TRUE)
expect_identical(is.matrix(out.IS$Z), TRUE)

expect_identical(dim(out.Frobenius$Z), dim(X))
expect_identical(dim(out.KL$Z), dim(X))
expect_identical(dim(out.IS$Z), dim(X))

### Test O-6: w

expect_identical(is.vector(out.Frobenius$w), TRUE)
expect_identical(is.vector(out.KL$w), TRUE)
expect_identical(is.vector(out.IS$w), TRUE)

expect_identical(length(out.Frobenius$w), 1)
expect_identical(length(out.KL$w), 1)
expect_identical(length(out.IS$w), 1)

### Test O-7: RecError

expect_identical(is.vector(out.Frobenius$RecError), TRUE)
expect_identical(is.vector(out.KL$RecError), TRUE)
expect_identical(is.vector(out.IS$RecError), TRUE)

### Test O-8 TrainRecError

expect_identical(is.vector(out.Frobenius$TrainRecError), TRUE)
expect_identical(is.vector(out.KL$TrainRecError), TRUE)
expect_identical(is.vector(out.IS$TrainRecError), TRUE)

### Test O-9: TestsRecError

expect_identical(is.vector(out.Frobenius$TestsRecError), TRUE)
expect_identical(is.vector(out.KL$TestsRecError), TRUE)
expect_identical(is.vector(out.IS$TestsRecError), TRUE)

### Test O-10: RelChange

expect_identical(is.vector(out.Frobenius$RelChange), TRUE)
expect_identical(is.vector(out.KL$RelChange), TRUE)
expect_identical(is.vector(out.IS$RelChange), TRUE)


## Test Error
### Test E-1: X

expect_error(ZITNMF(as.data.frame(X), J=J))

### Test E-2: Z

expect_error(ZITNMF(X, Z=cbind(X, X)))

### Test E-3: pseudocount

expect_error(ZITNMF(X, pseudocount="0.1"))
expect_error(ZITNMF(X, pseudocount=-0.1))

### Test E-4: initF

expect_error(ZITNMF(X, initF=X))
expect_error(ZITNMF(X, initF=NA))

### Test E-5: initA

expect_error(ZITNMF(X, initA=X))
expect_error(ZITNMF(X, initA=NA))

### Test E-6: fixF

expect_error(ZITNMF(X, fixF=X))
expect_error(ZITNMF(X, fixF=NA))

### Test E-7: fixA

expect_error(ZITNMF(X, fixA=X))
expect_error(ZITNMF(X, fixA=NA))

### Test E-8: init

expect_error(ZITNMF(X, init="NMFF"))

### Test E-9: J

expect_error(ZITNMF(X, J="5"))
expect_error(ZITNMF(X, J=c(2,4)))
expect_error(ZITNMF(X, J=10^10)

### Test E-10: Beta

expect_error(ZITNMF(X, J=J, Beta="0.1"))
expect_error(ZITNMF(X, J=J, Beta=TRUE))

### Test E-11: phi

expect_error(ZITNMF(X, J=J, phi="0.1"))

### Test E-12: thr

expect_error(ZITNMF(X, J=J, thr="0.1"))
expect_error(ZITNMF(X, J=J, thr=-2.3))

### Test E-13: num.iter

expect_error(ZITNMF(X, J=J, num.iter="100"))
expect_error(ZITNMF(X, J=J, num.iter=-1))

### Test E-14: verbose

expect_error(ZITNMF(X, J=J, verbose="verbose"))


## Test Decrease of Error

### Test D-1: RecError
.sampleRank <- function(x){
	rank(c(x[2], median(x), rev(x)[1]))
}

expect_identical(.sampleRank(out.Frobenius$RecError), 3:1)
expect_identical(.sampleRank(out.KL$RecError), 3:1)
expect_identical(.sampleRank(out.IS$RecError), 3:1)

### Test D-2: RelChange

expect_identical(.sampleRank(out.Frobenius$RelChange), 3:1)
expect_identical(.sampleRank(out.KL$RelChange), 3:1)
expect_identical(.sampleRank(out.IS$RelChange), 3:1)

## Session Information

sessionInfo()
