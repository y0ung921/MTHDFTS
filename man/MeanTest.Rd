% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MeanTest.R
\name{MeanTest}
\alias{MeanTest}
\title{Testing the mean functions for high-dimensional functional time series}
\usage{
MeanTest(p, M, Y, d = 3, prelist = seq(2, 11, length.out = 10), B)
}
\arguments{
\item{p}{The dimension of functional time series.}

\item{M}{The predetermined number of discrete points\eqn{(u_m(1\leq m\leq M))} in the functional domain\eqn{[0,1]}.}

\item{Y}{A \eqn{n*H}  matrix \eqn{(Y_{t,l})_{n*H}} with \eqn{H=p*M}; \eqn{Y_{t,l}=X_{t,j}(u_m)} with \eqn{j=\floor(l/M)} and  \eqn{m= l-(j-1)*M}.}

\item{d}{The number of principal components. The default value is 3.}

\item{prelist}{A set of the  blocksize. The default set is {2,3,4,...,11}.}

\item{B}{The number of bootstrap replications.}
}
\value{
An object of class "MeanTest" is a \eqn{1*2} matrix R,
where R[1,1] is the value of test statistic, R[1,2] is the corresponding p-value.
}
\description{
\code{MeanTest()} implements a new mean test proposed in
 Yang, Feng and Jiang (2024+) for the following hypothesis testing problem:
\deqn{H_0:\boldsymbol{\mu}(u)=0\mathrm{\ for\ any\ }u\in\mathcal{U},\mathrm{\ \ versus\ \ }
H_1:\boldsymbol{\mu}(u)\neq 0\mathrm{\ for\ some\ }u\in\mathcal{U}. }
}
\examples{
## Model
gendata = function(n,p,M){
S <- matrix(NA, nrow=p, ncol = p)
 for(i in 1:(p)){
   for(j in 1:(p)){
     S[i,j] = 0.6^(abs(i-j))
   }
 }
 eigenvalue <- eigen(S)$values
 eigenvector <- eigen(S)$vectors
 SS <- eigenvector \%*\% diag(sqrt(eigenvalue))\%*\% t(eigenvector)

 x <- NULL
 for (nn in 1:n) {
   xr.tmp <- matrix(rnorm(p*M,0,1), nrow = M, ncol = p)
   xr <- xr.tmp \%*\% SS
   BM.tmp <- apply(xr, 2, cumsum)/sqrt(M)
   x <- rbind(x, matrix(BM.tmp, nrow=1))
 }

 return(x)
}


library(fda)
## mean test
n <- 100
p <- 5
M <- 30
B <- 1000
Y <- gendata(n,p,M)
MeanTest(p,M,Y,d=3,prelist=seq(2,11,length.out =10),B )
}
\references{
Yang, L., Feng, Z. and Jiang, Q. (2024+).
\emph{Test for the mean of high-dimensional functional time series}.
}
