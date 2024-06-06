#' @name SmoothCurve
#' @title Local linear smooth for partially observed data
#' @description \code{SmoothCurve()} implements the local linear smoothing method.
#' @param p The dimension of functional time series.
#' @param M The predetermined number of discrete points\eqn{(u_m(1\leq m\leq M))} in the functional domain\eqn{[0,1]}.
#' @param Y a list of \eqn{n}, with each is a vector  \eqn{(X^*_{t,1,1},...,
#' X^*_{t,1,N_{t,1}},X^*_{t,2,1},\ldots,X^*_{t,2,N_{t,2}},\ldots,X^*_{t,p,N_{t,p}})}.
#' @param randp a list of \eqn{n}, with each is a vector  \eqn{(u^t_{1,1},...,
#' u^t_{1,N_{t,1}},u^t_{2,1},...,u^t_{2,N_{t,2}},...,u^t_{p,N_{t,p}})}.
#' @param N_vec a list of \eqn{n}, with each is a vector  \eqn{(N_{t,1},...,N_{t,p})}.
#'
#' @return An object of class "SmoothCurve" is a \eqn{n*H}  matrix
#' \eqn{(hatX_{t,v})_{n*H}} with \eqn{H=p*M}; \eqn{hatX_{t,v}=
#' hatX_{t,j}(u_m)} with \eqn{j=\floor(v/M)} and  \eqn{m= v-(j-1)*M}.
#'
#' @references Yang, L., Feng, Z. and Jiang, Q. (2024+).
#' \emph{Test for the mean of high-dimensional functional time series}.
#'
#' @examples
#'
#' # Assume N_{t,j}=N for all (t,j)
#' n <- 100
#' p <- 5
#' M <- 30
#' N <- 50
#' randp <- list()
#' for(t in 1:n){
#'   vec_test <- rep(NA,p*N)
#'   for(j in 1:p){
#'     vec_test[((j-1)*N+1):(j*N) ]<- sort(runif(N))
#'   }
#'   randp[[t]] <- vec_test
#' }
#' N_vec <- list()
#' for(t in 1:n){
#'   N_vec[[t]] <- rep(N,p)
#' }
#' Y <- list()
#' for(t in 1:n){
#'   Y[[t]] <- rnorm(p*N)
#' }
#' SmoothCurve(p,M,Y,randp,N_vec)
#'
#' @import fda
#' @import stats
#' @export
SmoothCurve <- function(p,M,Y,randp,N_vec){
  dis_p <- 1:M/M
  n <- length(Y)
  data_rec <- matrix(NA,nrow=n,ncol=p*M)
  for(t in 1:n){
    N_vec_new <- c(0,N_vec[[t]])
    for(j in 1:p){
      N_tj <- N_vec[[t]][j]
      h_tj <- N_tj^{-1/3}
      for(m in 1:M){
        aa <- exp(- ((randp[[t]][(sum(N_vec_new[1:j])+1): sum(N_vec_new[1:(j+1)])  ] - dis_p[m] )/h_tj)^2/2 )/(sqrt(2*pi)*0.6826895)
        bb1 <- (randp[[t]][(sum(N_vec_new[1:j])+1): sum(N_vec_new[1:(j+1)])  ] - dis_p[m] )/h_tj

        S_0 <- mean(aa)/h_tj
        S_1 <- mean(aa*bb1)/h_tj
        S_2 <- mean(aa*( randp[[t]][(sum(N_vec_new[1:j])+1): sum(N_vec_new[1:(j+1)])  ] - dis_p[m]  )^2/(h_tj^2))/h_tj

        R_0 <- mean(aa*Y[[t]][(sum(N_vec_new[1:j])+1): sum(N_vec_new[1:(j+1)])  ])/h_tj
        R_1 <- mean(aa*bb1*Y[[t]][(sum(N_vec_new[1:j])+1): sum(N_vec_new[1:(j+1)])  ])/h_tj

        data_rec[t,((j-1)*M)+m] <- (R_0*S_2-R_1*S_1)/(S_0*S_2-S_1^2+N_tj^{-2})
      }
    }
  }
  return(data_rec)
}



