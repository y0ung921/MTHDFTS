#' @name MeanTest
#' @title Testing the mean functions for high-dimensional functional time series
#' @description \code{MeanTest()} implements a new mean test proposed in
#'  Yang, Feng and Jiang (2024+) for the following hypothesis testing problem:
#' \deqn{H_0:\boldsymbol{\mu}(u)=0\mathrm{\ for\ any\ }u\in\mathcal{U},\mathrm{\ \ versus\ \ }
#' H_1:\boldsymbol{\mu}(u)\neq 0\mathrm{\ for\ some\ }u\in\mathcal{U}. }
#' @param p The dimension of functional time series.
#' @param M The predetermined number of discrete points\eqn{(u_m(1\leq m\leq M))} in the functional domain\eqn{[0,1]}.
#' @param Y A \eqn{n*H}  matrix \eqn{(Y_{t,l})_{n*H}} with \eqn{H=p*M}; \eqn{Y_{t,l}=X_{t,j}(u_m)} with \eqn{j=\floor(l/M)} and  \eqn{m= l-(j-1)*M}.
#' @param d The number of principal components. The default value is 3.
#' @param prelist A set of the  blocksize. The default set is {2,3,4,...,11}.
#' @param B The number of bootstrap replications.
#' @return An object of class "MeanTest" is a \eqn{1*2} matrix R,
#' where R[1,1] is the value of test statistic, R[1,2] is the corresponding p-value.
#' @references Yang, L., Feng, Z. and Jiang, Q. (2024+).
#' \emph{Test for the mean of high-dimensional functional time series}.
#' @examples
#' ## Model
#' gendata = function(n,p,M){
#' S <- matrix(NA, nrow=p, ncol = p)
#'  for(i in 1:(p)){
#'    for(j in 1:(p)){
#'      S[i,j] = 0.6^(abs(i-j))
#'    }
#'  }
#'  eigenvalue <- eigen(S)$values
#'  eigenvector <- eigen(S)$vectors
#'  SS <- eigenvector %*% diag(sqrt(eigenvalue))%*% t(eigenvector)
#'
#'  x <- NULL
#'  for (nn in 1:n) {
#'    xr.tmp <- matrix(rnorm(p*M,0,1), nrow = M, ncol = p)
#'    xr <- xr.tmp %*% SS
#'    BM.tmp <- apply(xr, 2, cumsum)/sqrt(M)
#'    x <- rbind(x, matrix(BM.tmp, nrow=1))
#'  }
#'
#'  return(x)
#'}
#'
#'
#' library(fda)
#' ## mean test
#' n <- 100
#' p <- 5
#' M <- 30
#' B <- 1000
#' Y <- gendata(n,p,M)
#' MeanTest(p,M,Y,d=3,prelist=seq(2,11,length.out =10),B )
#' @import fda
#' @import stats
#' @export
MeanTest <- function(p,M,Y,d=3,prelist=seq(2,11,length.out =10),B ){
  ## Projection
  n <- dim(Y)[1]
  Y_cen <- t(t(Y) - colMeans(Y))
  basisobj = create.fourier.basis(c(0,1), nbasis=49)
  basismatrix = eval.basis(1:M/M, basisobj)
  Y_pro <- matrix(NA, nrow=n, ncol= p*d)
  for(j in 1:p){
    dnm <- (t(Y_cen[,((j-1)*M+1):(j*M)])%*%Y_cen[,((j-1)*M+1):(j*M)])/n
    D1 <- (t(basismatrix)%*% dnm %*% basismatrix)/(M^2)
    d.eigen <- eigen(D1)
    evec <- Re(d.eigen$vectors[,1:d])
    for(i in 1:d){
      ree <- 0
      for(l in 1:49){
        ree <- ree+evec[l,i]*(  (Y[,((j-1)*M+1):(j*M)] %*% basismatrix[,l]))/M
      }
      Y_pro[,(j-1)*d+i] <- ree
    }
  }
  # The value of test statistic
  Tn <- max((colSums( Y_pro))^2 %*% (kronecker( diag(p), rep(1,d))))/n

  Y_procen <- t(t(Y_pro) - colMeans(Y_pro))

  blength <- length(prelist)
  blist <- c( 2*prelist[1]-prelist[2], prelist, 2*prelist[blength]-prelist[blength-1] )
  bli <- rep(NA,blength)
  for(s in 1:blength){
    ml1 <- blist[s]
    S1 <- ceiling(n/ml1)
    namatrix1 <- colSums((t(((kronecker(diag(rep(1,S1)),rep(1,ml1)))[1:n,]))%*%Y_procen)^2)
    ml2 <- blist[s+1]
    S2 <- ceiling(n/ml2)
    namatrix2 <- colSums((t(((kronecker(diag(rep(1,S2)),rep(1,ml2)))[1:n,]))%*%Y_procen)^2)
    ml3 <- blist[s+2]
    S3 <- ceiling(n/ml3)
    namatrix3 <- colSums((t(((kronecker(diag(rep(1,S3)),rep(1,ml3)))[1:n,]))%*%Y_procen)^2)

    ml <- rbind(namatrix1,namatrix2,namatrix3)
    bli[s] <- sum(apply(ml,2,sd))
  }
  lenopt <- prelist[which.min(bli)]

  S <- ceiling(n/lenopt)
  normal <- (kronecker(diag(rep(1,S)),rep(1,lenopt))%*%matrix(rnorm(B*S),nrow=S,ncol=B))[1:n,]

  Tn_G <- rep(NA,B)
  for(b in 1:B){
    comb <- Y_procen * normal[,b]
    Tn_G[b] <- max((colSums( comb))^2 %*% (kronecker( diag(p), rep(1,d))))/n
  }
  sortG <- sort(Tn_G)
  loc <- which(sortG>=Tn)[1]
  if(is.na(loc)==TRUE){
    p_value <- 1
  }else{
    p_value <- (B-loc)/B
  }
  re_matrix <- matrix(NA,nrow=1,ncol=2)
  colnames(re_matrix) <- c('T_n','p_value')
  re_matrix[1,1] <- Tn
  re_matrix[1,2] <- p_value
  return(re_matrix)
}
