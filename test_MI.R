# Test code of MI Function
## Dependent Packages Installation
if (!requireNamespace("testthat", quietly = TRUE)){
    install.packages("testthat", repos="https://cran.r-project.org")
}
library("testthat")
if (!requireNamespace("Rcpp", quietly = TRUE)){
    install.packages("Rcpp", repos="https://cran.r-project.org")
}
library("Rcpp")
if (!requireNamespace("infotheo", quietly = TRUE)){
    install.packages("infotheo", repos="https://cran.r-project.org")
}
library("infotheo")

## MI Function (Paste your MI code here)
MI <- function(X, pseudocount=1e-10,
    algorithm = c("fastpairmi", "hist", "knn"),
    h=0.1, # Window Size of Gauss Kernel used in "fastpairmi"
    size=20, # Number of histgrams used in "hist"
    k=5, # Number of neighborhoods used in "knn"
    # 必要に応じて適宜パラメーターを追加
	verbose=FALSE){
    ######################################
    # Argument Check
    ######################################
    .checkMI(X, h, size, k, verbose)
    algorithm <- match.arg(algorithm)
    ######################################
    # Initialization
    #######################################
    int <- .initMI(X)
    X <- int$X
    ######################################
    # MI Calculation
    #######################################
    .MI_METHODS[algorithm](X, h, size, k)
}

# 必要に応じて適宜パラメーターチェックを追加
.checkMI <- function(X, h, size, k, verbose){
    stopifnot(is.matrix(X))
    stopifnot(is.numeric(h))
    stopifnot(h > 0)
    stopifnot(is.numeric(size))
    stopifnot(size > 0)
    stopifnot(size %% 1 == 0)
    stopifnot(is.numeric(k))
    stopifnot(k > 0)
    stopifnot(k %% 1 == 0)
    stopifnot(is.logical(verbose))
}

.initMI <- function(X){
    # Xに前処理などあれば
}



.FastPairMI <- function(X, h, size, k){
    # FastPairMIの具体的な実装（C++コード呼び出し）
}

.HistMI <- function(X, h, size, k){
    # ヒストグラムベースのMIの具体的な実装
}

.kNNMI <- function(X, h, size, k){
    # kNNベースのMIの具体的な実装
}

.MI_METHODS <- list(
    "fastpairmi" = .FastPairMI,
    "hist" = .HistMI,
    "knn" = .kNNMI)

## Simulation Dataset
N <- 1000 # Number of Data (e.g., Gene)
M <- 300 # Number of Dimension (e.g., Cell Pairs)
X <- matrix(runif(N*M), nrow=N, ncol=M)

## Perform MI against Simulation Dataset
outFastPairMI <- MI(X, algorithm="fastpairmi")
outHistMI <- MI(X, algorithm="hist")
outkNNMI <- MI(X, algorithm="knn")

## Test Input object / type
### Test I-1: Object Names
expect_identical(names(formals(MI)),
    c("X", "algorithm", "h", "size", "k", "verbose"))

### Test I-2: X
expect_identical(as.character(formals(MI)$X), "")

### Test I-3: algorithm
expect_identical(formals(MI)$algorithm, c("fastpairmi", "hist", "knn"))

### Test I-4: h
expect_identical(formals(MI)$h, 0.1)

### Test I-5: size
expect_identical(formals(MI)$size, 20)

### Test I-6: k
expect_identical(formals(MI)$k, 5)

### Test I-7: verbose
expect_identical(formals(MI)$verbose, FALSE)

## Test Output object / type
### Test O-1: Object
expect_identical(is.list(outFastPairMI), TRUE)
expect_identical(is.list(outHistMI), TRUE)
expect_identical(is.list(outKNNMI), TRUE)

### Test O-2: Object Size
expect_identical(dim(outFastPairMI), c(N, N))
expect_identical(dim(outHistMI), c(N, N))
expect_identical(dim(outKNNMI), c(N, N))

## Test Error
### Test E-1: X
expect_error(MI(as.data.frame(X)))

### Test E-2: h
expect_error(MI(X, h="0.4"))
expect_error(MI(X, h=-2.3))

### Test E-2: size
expect_error(MI(X, size="20"))
expect_error(MI(X, size=-5))

### Test E-2: k
expect_error(MI(X, k="2"))
expect_error(MI(X, k=-3))

### Test E-2: verbose
expect_error(MI(X, verbose="verbose"))

## Test Speed
# 連続量データ
X_small <- matrix(runif(10*30*10), nrow=10, ncol=30*10)
X_medium <- matrix(runif(100*30*10), nrow=100, ncol=30*10)
X_large <- matrix(runif(1000*30*10), nrow=1000, ncol=30*10)

# 離散化済みデータを想定
dX_small <- matrix(rbinom(10*30*10,2,0.5), nrow=10, ncol=30*10)
dX_medium <- matrix(rbinom(100*30*10,2,0.5), nrow=100, ncol=30*10)
dX_large <- matrix(rbinom(1000*30*10,2,0.5), nrow=1000, ncol=30*10)

# MI by infortheo
dPairMI <- function(X){
	apply(X, 1, function(x) apply(X, 1, function(xx) mutinformation(x, xx)))
}

out_small <- system.time(dPairMI(dX_small))[3] # 0.028 (s)
out_medium <- system.time(dPairMI(dX_medium))[3] # 2.689 (s)
out_large <- system.time(dPairMI(dX_large))[3] # 270.726 (s)

# Test MI Algorithms' speed
expect_true(system.time(MI(X_small, algorithm="fastpairmi")) < out_small / 10)
expect_true(system.time(MI(X_medium, algorithm="fastpairmi")) < out_medium / 10)
expect_true(system.time(MI(X_large, algorithm="fastpairmi")) < out_large / 10)
expect_true(system.time(MI(X_small, algorithm="hist")) < out_small / 10)
expect_true(system.time(MI(X_medium, algorithm="hist")) < out_medium / 10)
expect_true(system.time(MI(X_large, algorithm="hist")) < out_large / 10)
expect_true(system.time(MI(X_small, algorithm="knn")) < out_small / 10)
expect_true(system.time(MI(X_medium, algorithm="knn")) < out_medium / 10)
expect_true(system.time(MI(X_large, algorithm="knn")) < out_large / 10)

## Test Accuracy
# ・でたらめな結果を出力しているわけではないことを示すテスト
# ・例: ランダムなデータの中に、非線形な依存関係がある変数ペアを埋め込んで、MIの値が実際に大きくなるか
# ・例: 原著論文で行われた検証の再現

## Session Information
sessionInfo()
