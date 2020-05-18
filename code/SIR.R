
SIR <- function (X, y, k, h = 10, natH=FALSE,...) 
{
    X <- as.matrix(X)
    n <- nrow(X)
    p <- ncol(X)
    MEAN <- colMeans(X)
    Xc <- sweep(X, 2, MEAN, "-")
    COV <- crossprod(Xc)/n
    EVD.COV <- eigen(COV, symmetric = TRUE)
    COV.inv.sqrt <- EVD.COV$vectors %*% diag((1/EVD.COV$values)^0.5) %*% 
        t(EVD.COV$vectors)
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