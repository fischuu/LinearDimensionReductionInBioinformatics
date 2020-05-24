SIR.natH <- function (X, y, k, h = 10, natH=FALSE,...) 
{
    X <- as.matrix(X)
    n <- nrow(X)
    p <- ncol(X)
    MEAN <- colMeans(X)
    Xc <- sweep(X, 2, MEAN, "-")
    COV <- crossprod(Xc)/n
    
    EVD.COV <- eigen(COV, symmetric = TRUE)
    COV.inv.sqrt <- EVD.COV$vectors %*% diag((1/EVD.COV$values)^0.5) %*% t(EVD.COV$vectors)

    Y <- tcrossprod(Xc, COV.inv.sqrt)
    if(natH){
      S2 <- covSIRnat(Y, y, ...)
    } else {
      S2 <- covSIR(Y, y, h = h, ...)
    }
    EVD.S2 <- eigen(S2, symmetric = TRUE)
    D <- EVD.S2$values
    W <- crossprod(EVD.S2$vectors, COV.inv.sqrt)
    S <- tcrossprod(Xc, W)
    colnames(S) <- paste0("SIC.", 1:p)
    RES <- list(W = W, S = S, D = D, MU = MEAN)
    RES
}

SIRasymp.natH <- function(X, y, k, h=10, natH=FALSE, ...)
{
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  MEAN <- colMeans(X)
  Xc <- sweep(X, 2, MEAN, "-")
  COV <- crossprod(Xc)/n
  
  EVD.COV <- eigen(COV, symmetric=TRUE)
  COV.inv.sqrt <- EVD.COV$vectors %*% diag((1/EVD.COV$values)^0.5) %*%  t(EVD.COV$vectors)
  
  Y <- tcrossprod(Xc, COV.inv.sqrt)
  if(natH){
    S2 <- covSIRnat(Y, y, ...)
  } else {
    S2 <- covSIR(Y, y, h = h, ...)
  }
  EVD.S2 <- eigen(S2, symmetric = TRUE)
  
  D <- EVD.S2$values
  W <- crossprod(EVD.S2$vectors, COV.inv.sqrt)
  S <- tcrossprod(Xc, W)
  
  TESTSTATISTIC <-  n*sum(D[(k+1):p])
  names(TESTSTATISTIC) = "T"
  
  PARAMETER <- (p-k)*(h-k-1)
  names(PARAMETER) <- "df"
  PVAL <- 1- pchisq(TESTSTATISTIC, df=PARAMETER)
  METHOD <- "SIR test for subspace dimension"
  ALTERNATIVE <- paste0("the last ",p-k, " eigenvalues are not zero")
  
  colnames(S) <- paste0("SIC.",1:p)
  
  RES <- list(statistic = TESTSTATISTIC, p.value = PVAL, parameter = PARAMETER, method=METHOD, alternative = ALTERNATIVE, k=k, W=W, S=S, D=D,
              MU=MEAN)
  class(RES) <- c("ictest", "htest")
  RES
}

# internal:
# Function to compute the test statistic

SIR_boot_teststatistic <- function(X, y, k, h, natH, ...)
{
  n <- nrow(X)
  p <- ncol(X)
  MEAN <- colMeans(X)
  Xc <- sweep(X, 2, MEAN, "-")
  COV <- crossprod(Xc)/n
  
  EVD.COV <- eigen(COV, symmetric=TRUE)
  COV.inv.sqrt <- EVD.COV$vectors %*% diag((1/EVD.COV$values)^0.5) %*%  t(EVD.COV$vectors)
  
  Y <- tcrossprod(Xc, COV.inv.sqrt)
  if(natH){
    S2 <- covSIRnat(Y, y, ...)
  } else {
    S2 <- covSIR(Y, y, h = h, ...)
  }
  EVD.S2 <- eigen(S2, symmetric = TRUE)
  
  D <- EVD.S2$values
  
  W <- crossprod(EVD.S2$vectors, COV.inv.sqrt)
  
  TEST.STATISTIC.X <- sum(D[(k+1):p])
  return(TEST.STATISTIC.X)
}

# internal:
# Function to create the bootstrap sample and to compute its test statistic

SIR_boot_normal_all <- function(Z1,Z2,y, Winv, MEAN, k, h, ...)
{
  n <- nrow(Z1)
  
  ind <- sample(1:n, replace=TRUE)
  
  Z1tilde <- Z1[ind,]
  Z2tilde<-Z2[sample(1:n,replace=TRUE),]
  Zstar <- cbind(Z1tilde,Z2tilde)
  
  Xstar <- tcrossprod(Zstar, Winv)
  Xstar <- sweep(Xstar, 2, MEAN, "+")
  
  ystar <- y[ind]
  TEST.STATISTIC.Xstar <- SIR_boot_teststatistic(Xstar, ystar, k, h=h, ...)
  TEST.STATISTIC.Xstar
}

# General function for SIR bootstrapping:

SIRboot.natH <- function(X, y, k, h=10, n.boot = 200, natH=FALSE,...)
{
  DNAME <- deparse(substitute(X))
  n <- nrow(X)
  p <- ncol(X)
  MEAN <- colMeans(X)
  Xc <- sweep(X, 2, MEAN, "-")
  COV <- crossprod(Xc)/n
  
  EVD.COV <- eigen(COV, symmetric=TRUE)
  COV.inv.sqrt <- EVD.COV$vectors %*% diag((1/EVD.COV$values)^0.5) %*%  t(EVD.COV$vectors)
  
  Y <- tcrossprod(Xc, COV.inv.sqrt)
  if(natH){
    S2 <- covSIRnat(Y, y, ...)
  } else {
    S2 <- covSIR(Y, y, h = h, ...)
  }
  EVD.S2 <- eigen(S2, symmetric = TRUE)
  
  D <- EVD.S2$values
  W <- crossprod(EVD.S2$vectors, COV.inv.sqrt)
  
  TEST.STATISTIC.X <-  sum(D[(k+1):p])
  
  Z <- tcrossprod(Xc, W)
  
  Z1 <- Z[,0:k, drop=FALSE]
  Z2 <- Z[,(k+1):p, drop=FALSE]
  
  Winv <- solve(W)
  
  TEST.STATISTICS.Xstar <- replicate(n.boot, SIR_boot_normal_all(Z1, Z2, y,  Winv, MEAN, k, h, ...))
  PVAL <- (sum(TEST.STATISTIC.X<TEST.STATISTICS.Xstar)+1)/(n.boot+1)
  METHOD <- "SIR bootstrapping test for subspace dimension"
  ALTERNATIVE <- paste0("the last ",p-k, " eigenvalues are not zero")
  PARAMETER <- n.boot
  names(PARAMETER) <-  "replications"
  names(TEST.STATISTIC.X) <- "T"
  
  colnames(Z) <- paste0("SIC.",1:p)
  
  RES <- list(statistic = TEST.STATISTIC.X, p.value = PVAL, parameter = PARAMETER, method=METHOD, data.name=DNAME, alternative = ALTERNATIVE, k=k, W=W, S=Z, D=D,
              MU=MEAN)
  class(RES) <- c("ictest", "htest")
  RES
}



 covSIR <- function (X, y, h = 10, ...) 
{
    n <- nrow(X)
    slices <- as.matrix(cut(y, breaks = c(quantile(y, probs = seq(0, 
        1, by = 1/h), ...)), include.lowest = TRUE, labels = FALSE))
    slicemean <- do.call(rbind, by(X, slices, colMeans))
    cov(slicemean)
 }
 
 covSIRnat <- function (X, y, ...) 
 {
   n <- nrow(X)
   slices <- as.matrix(y)
   slicemean <- do.call(rbind, by(X, slices, colMeans))
   cov(slicemean)
 }