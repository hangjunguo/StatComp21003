#' @title Simluate Toeplitz design and its response
#' @description Simluate Toeplitz design and its response
#' @param n the datasize
#' @param p the number of variables
#' @param rho the correlation between variables
#' @param truebeta the regression coefficients
#' @return Toeplitz design and its response
#' @examples
#' \dontrun{
#' n <- 200;  p <- 300;  rho <- 0.8; s <- 20
#' S_0 <- c(1:s)
#' truebeta <- rep(0, p)
#' truebeta[S_0] <- runif(s, 0, 2)
#' dat <- simdata(n, p, rho, truebeta)
#' }
#' @export
simdata <- function(n, p, rho, truebeta){
  
  set.seed(99)
  covmat <- diag(p)
  for(i in 1:(p-1)){
    for(j in (i+1):p){
      covmat[i, j] <- rho^abs(i-j)
    }
  }
  covmat <- covmat+t(covmat)-diag(p)
  x_design <- mvrnorm(n, mu = rep(0, p), Sigma = covmat)
  y_response <- x_design %*% truebeta + rnorm(n)
  
  return(list(X=x_design, y=y_response))
}



#' @title Low Dimensional Projection Estimation
#' @description Algorithms for fitting Low Dimensional Projection Estimation
#' @param X the design matrix
#' @param y the response
#' @param alpha the significance level alpha, default is 0.05
#' @return A LDPE object is returned, including confidence intervals for regression coefficients. 
#' @examples
#' \dontrun{
#' dat <- simdata(n, p, rho, truebeta)
#' X <- dat$X
#' y <- dat$y
#' ldpe <- LDP(X, y)
#' CI <- ldpe$CI; len_CI <- ldpe$len_CI
#' }
#' @import microbenchmark
#' @import Rcpp
#' @import stats
#' @import ggplot2
#' @import MASS
#' @import bootstrap
#' @import boot
#' @import glmnet
#' @import ncpen
#' @import SIS
#' @import VariableScreening
#' @import bestglm
#' @import hdi
#' @import scalreg
#' @useDynLib WAP
#' @export
LDP <- function(X, y, alpha=0.05){
  
  n <- nrow(X); p <- ncol(X)
  CI <- matrix(nrow = 2, ncol = p)
  len_CI <- numeric(p)
  z <- matrix(nrow = n, ncol = p)
  
  ###initial scaled lasso estimator
  scallasso <- scalreg(X, y, lam0 = "univ")
  hsigma <- scallasso$hsigma
  beta <- as.matrix(scallasso$coefficients)
  
  for(j in 1:p){
    
    ###nodewise lasso, tuned by 10-fold CV
    lasso1 <- glmnet(X[, -j], X[, j], family = 'gaussian',
                     intercept = FALSE, standardize = TRUE,
                     alpha = 1, nlambda = 100)
    omega <- as.numeric(coef(lasso1, s=lasso1$lambda[100]))[-1]
    e_lasso1 <- X[, j] - X[, -j]%*%omega
    
    ###relaxed orthogonal vector z_j
    z[ ,j] <- e_lasso1
    debias <- as.numeric((t(z[ ,j])%*%(y-X%*%beta)) / (t(z[ ,j])%*%X[, j]))
    beta[j] <- beta[j] + debias
    tau <- as.numeric(sqrt(t(z[ ,j]%*%z[ ,j])) / abs(t(z[ ,j])%*%X[, j]))
    len_CI[j] <- 2*qnorm(1-alpha/2)*tau*hsigma
    CI[, j] <- c(beta[j]-len_CI[j]/2, beta[j]+len_CI[j]/2)
    
  }
  
  return(list(z=z, beta=beta, CI=CI, len_CI=len_CI))
  
}



#' @title Weighted Adaptive Projection
#' @description Weighted Adaptive Projection
#' @param X the design matrix
#' @param y the response
#' @param S the pre-screened identifiable feature set
#' @param alpha the significance level alpha, default is 0.05
#' @return A WAP object is returned, including the estimated confidence intervals for LR model parameters. 
#' @examples
#' \dontrun{
#' dat <- simdata(n, p, rho, truebeta)
#' X <- dat$X
#' y <- dat$y
#' S <- SIS(X, y, family = 'gaussian', tune = 'bic', seed = 2)$sis.ix0
#' wap <- WAP(X, y, S)
#' CI <- wap$CI; len_CI <- wap$len_CI
#' }
#' @export
WAP <- function(X, y, S, alpha=0.05){
  #X <- x_design; y <- y_response; 
  #alpha=0.05; set.seed(99)
  n <- nrow(X); p <- ncol(X)
  CI <- matrix(nrow = 2, ncol = p)
  len_CI <- numeric(p)
  z <- matrix(nrow = n, ncol = p)
  
  
  ###initial scaled lasso estimator
  scallasso <- scalreg(X, y, lam0 = "univ")
  hsigma <- scallasso$hsigma
  beta <- scallasso$coefficients
  
  #delicate isolation
  DeliIsolate <- function(j, S, X){
    
    n <- nrow(X); p <- ncol(X)
    
    stdnorm <- function(x){
      x_stdnorm <- sqrt(sum(x^2)) / sqrt(length(x))
      return(x_stdnorm)
    }
    
    if(j %in% S){
      jj <- which(S==j)
      Pmat <- X[,S[-jj]] %*% solve(t(X[,S[-jj]])%*%X[,S[-jj]]) %*% t(X[,S[-jj]])
      S_cj <- c(1:p)[-S[-jj]]
      
      psi <- matrix(nrow=n, ncol=length(S_cj))
      v <- rep(0, length(S_cj))
      counter <- 1
      for(k in S_cj){
        psi[, counter] <- (diag(n)-Pmat) %*% X[,k]
        v[counter] <- stdnorm(psi[, counter])
        counter <- counter+1
      }
    }else{
      S_c <- c(1:p)[-S]
      Pmat <- X[,S] %*% solve(t(X[,S])%*%X[,S]) %*% t(X[,S])
      psi <- (diag(n)-Pmat) %*% X[,S_c]
      v <- apply(psi, 2, stdnorm)
    }
    
    return(list(psi=psi, v=v))
  }
  
  #relaxed projection
  RelxProj <- function(j, S, psi, v, beta, alpha=0.05){
    ###find the location of the target tg
    if(j %in% S){
      jj <- which(S==j)
      S_cj <- c(1:p)[-S[-jj]]
      tg <- which(S_cj==j)
    }else{
      S_c <- c(1:p)[-S]
      tg <- which(S_c==j)
    }
    
    
    ###nodewise lasso, tuned by GIC
    lasso2 <- ncpen(psi[, tg], psi[, -tg],
                    family = 'gaussian', penalty = 'lasso',
                    intercept = FALSE, x.standardize = TRUE)
    omega <- gic.ncpen(lasso2, verbose = FALSE)$opt.beta
    e_lasso2 <- psi[, tg] - psi[, -tg]%*%omega
    
    ###adaptive orthogonal vector z_j
    z <- e_lasso2
    beta <- as.numeric((t(z)%*%y) / (t(z)%*%X[, j]))
    tau <- as.numeric(sqrt(t(z)%*%z) / abs(t(z)%*%X[, j]))
    len_CI <- 2*qnorm(1-alpha/2)*tau*hsigma
    CI <- c(beta-len_CI/2, beta+len_CI/2)
    
    
    return(list(z=z, beta=beta, CI=CI, len_CI=len_CI))
  }
  
  
  
  for(j in S){
    
    psi <- DeliIsolate(j, S, X)$psi
    v <- DeliIsolate(j, S, X)$v
    
    rp <- RelxProj(j, S, psi, v, beta[j])
    z[, j] <- rp$z
    beta[j] <- rp$beta
    CI[, j] <- rp$CI
    len_CI[j] <- rp$len_CI
    
  }
  
  psi <- DeliIsolate(-1, S, X)$psi#psi is same for j in S_c
  v <- DeliIsolate(-1, S, X)$v
  S_c <- c(1:p)[-S]
  
  for(j in S_c){
    
    rp <- RelxProj(j, S, psi,v, beta[j])
    z[, j] <- rp$z
    beta[j] <- rp$beta
    CI[, j] <- rp$CI
    len_CI[j] <- rp$len_CI
    
  }
  
  return(list(CI=CI, len_CI=len_CI, beta=beta, z=z))
}

