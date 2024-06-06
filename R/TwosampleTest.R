#' @name TwosampleTest
#' @title Two sample testing function for high-dimensional functional time series
#' @description \code{TwosampleTest()} implements a new two-sample test proposed in
#'  Yang, Feng and Jiang (2024+) for the following hypothesis testing problem:
#' \deqn{H_0:\boldsymbol{\mu}_1(u)=\boldsymbol{\mu}_2(u)\mathrm{\ for\ any\ }u\in\mathcal{U},\mathrm{\ \ versus\ \ }
#' H_1:\boldsymbol{\mu}_1(u)\neq \boldsymbol{\mu}_2(u)\mathrm{\ for\ some\ }u\in\mathcal{U}. }
#' @param p The dimension of functional time series.
#' @param M The predetermined number of discrete points\eqn{(u_m(1\leq m\leq M))} in the functional domain\eqn{[0,1]}.
#' @param Y1 A \eqn{n*H}  matrix \eqn{(Y1_{t,l})_{n*H}} with \eqn{H=p*M}; \eqn{Y1_{t,l}=X1_{t,j}(u_m)} with \eqn{j=\floor(l/M)} and  \eqn{m= l-(j-1)*M}.
#' @param Y2 A \eqn{n*H}  matrix \eqn{(Y2_{t,l})_{n*H}} with \eqn{H=p*M}; \eqn{Y2_{t,l}=X2_{t,j}(u_m)} with \eqn{j=\floor(l/M)} and  \eqn{m= l-(j-1)*M}.
#' @param d The number of principal components. The default value is 3.
#' @param prelist A set of the  blocksize. The default set is {2,3,4,...,11}.
#' @param B The number of bootstrap replications.
#'
#' @return An object of class "TwosampleTest" is a \eqn{1*2} matrix R,
#' where R[1,1] is the value of test statistic, R[1,2] is the corresponding p-value.
#'
#' @references Yang, L., Feng, Z. and Jiang, Q. (2024+).
#' \emph{Test for the mean of high-dimensional functional time series}.
#'
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
#' ## two sample test
#' n_1 <- 100
#' n_2 <- 120
#' p <- 5
#' M <- 30
#' B <- 1000
#' Y1 <- gendata(n_1,p,M)
#' Y2 <- gendata(n_2,p,M)
#' TwosampleTest(p,M,Y1,Y2,d=3,prelist=seq(2,11,length.out =10),B)
#' @import fda
#' @import stats

#' @export
TwosampleTest <- function(p,M,Y1,Y2,d=3,prelist=seq(2,11,length.out =10),B){
  ## Projection
  n_1 <- dim(Y1)[1]
  n_2 <- dim(Y2)[1]
  theta <- n_1/(n_1+n_2)
  Y1_cen <- t(t(Y1) - colMeans(Y1))
  Y2_cen <- t(t(Y2) - colMeans(Y2))
  basisobj = create.fourier.basis(c(0,1), nbasis=49)
  basismatrix = eval.basis(1:M/M, basisobj)
  Y1_pro <- matrix(NA, nrow=n_1, ncol= p*d)
  Y2_pro <- matrix(NA, nrow=n_2, ncol= p*d)
  for(j in 1:p){
    CN <- (t(Y1_cen[,((j-1)*M+1):(j*M)])%*%Y1_cen[,((j-1)*M+1):(j*M)])/n_1
    CM <- (t(Y2_cen[,((j-1)*M+1):(j*M)])%*%Y2_cen[,((j-1)*M+1):(j*M)])/n_2
    dnm <- (1-theta)*CN+theta*CM
    D <- (t(basismatrix)%*% dnm %*% basismatrix)/(M^2)
    d.eigen <- eigen(D)
    evec <- Re(d.eigen$vectors[,1:d])
    for(i in 1:d){
      ree1 <- 0
      ree2 <- 0
      for(l in 1:49){
        ree1 <- ree1+evec[l,i]*(  (Y1[,((j-1)*M+1):(j*M)] %*% basismatrix[,l]))/M
        ree2 <- ree2+evec[l,i]*(  (Y2[,((j-1)*M+1):(j*M)] %*% basismatrix[,l]))/M
      }
      Y1_pro[,(j-1)*d+i] <- ree1
      Y2_pro[,(j-1)*d+i] <- ree2
    }
  }
  # The value of test statistic
  Tn_ts <- max((colMeans(Y1_pro)-colMeans(Y2_pro))^2%*% (kronecker( diag(p), rep(1,d))))*(n_1*n_2/(n_1+n_2))


  Y1_cen <- t(t(Y1_pro) - colMeans(Y1_pro))
  Y2_cen <- t(t(Y2_pro) - colMeans(Y2_pro))

  blength <- length(prelist)
  blist <- c( 2*prelist[1]-prelist[2], prelist, 2*prelist[blength]-prelist[blength-1] )
  bli <- matrix(NA,nrow=blength,ncol=blength)
  for(s1 in 1:blength){
    for(s2 in 1:blength){
      ml11 <- blist[s1]
      S11 <- ceiling(n_1/ml11)
      ml21 <- blist[s1+1]
      S21 <- ceiling(n_1/ml21)
      ml31 <- blist[s1+2]
      S31 <- ceiling(n_1/ml31)
      ml12 <- blist[s2]
      S12 <- ceiling(n_2/ml12)
      ml22 <- blist[s2+1]
      S22 <- ceiling(n_2/ml22)
      ml32 <- blist[s2+2]
      S32 <- ceiling(n_2/ml32)
      namatrix1 <- colSums((t(((kronecker(diag(rep(1,S11)),rep(1,ml11)))[1:n_1,]))%*%Y1_cen)^2)*(1-theta)/n_1+colSums((t(((kronecker(diag(rep(1,S12)),rep(1,ml12)))[1:n_2,]))%*%Y2_cen)^2)*(theta)/n_2
      namatrix2 <- colSums((t(((kronecker(diag(rep(1,S21)),rep(1,ml21)))[1:n_1,]))%*%Y1_cen)^2)*(1-theta)/n_1+colSums((t(((kronecker(diag(rep(1,S22)),rep(1,ml22)))[1:n_2,]))%*%Y2_cen)^2)*(theta)/n_2
      namatrix3 <- colSums((t(((kronecker(diag(rep(1,S31)),rep(1,ml31)))[1:n_1,]))%*%Y1_cen)^2)*(1-theta)/n_1+colSums((t(((kronecker(diag(rep(1,S32)),rep(1,ml32)))[1:n_2,]))%*%Y2_cen)^2)*(theta)/n_2

      ml <- rbind(namatrix1,namatrix2,namatrix3)
      bli[s1,s2] <- sum(apply(ml,2,sd))
    }
  }
  lenopt1 <- as.numeric(which(bli==min(bli),arr.ind=TRUE)[1,1])
  lenopt2 <- as.numeric(which(bli==min(bli),arr.ind=TRUE)[1,2])

  S1 <- ceiling(n_1/lenopt1)
  S2 <- ceiling(n_2/lenopt2)
  normal1 <- (kronecker(diag(rep(1,S1)),rep(1,lenopt1))%*%matrix(rnorm(B*S1),nrow=S1,ncol=B))[1:n_1,]
  normal2 <- (kronecker(diag(rep(1,S2)),rep(1,lenopt2))%*%matrix(rnorm(B*S2),nrow=S2,ncol=B))[1:n_2,]

  Tn_ts_G <- rep(NA,B)
  for(b in 1:B){
    comb1 <- Y1_cen * normal1[,b]
    comb2 <- Y2_cen * normal2[,b]
    Tn_ts_G[b] <- max((colMeans( comb1)-colMeans( comb2))^2 %*% (kronecker( diag(p), rep(1,d))))*(n_1*n_2)/(n_1+n_2)
  }
  sortG_ts <- sort(Tn_ts_G)
  loc_ts <- which(Tn_ts_G>=Tn_ts)[1]
  if(is.na(loc_ts)==TRUE){
    p_value <- 1
  }else{
    p_value <- (B-loc_ts)/B
  }
  re_matrix <- matrix(NA,nrow=1,ncol=2)
  colnames(re_matrix) <- c('T_n','p_value')
  re_matrix[1,1] <- Tn_ts
  re_matrix[1,2] <- p_value
  return(re_matrix)
}


