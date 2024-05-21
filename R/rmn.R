rmn <- function(k, M, U, V) {

  n <- dim(M)[1]  ;  p <- dim(M)[2]
  A <- chol(U)  ;  B <- chol(V)
  Y <- list()
  for ( i in 1:k ) {
    X <- Rfast::matrnorm(n, p)
    Y[[ i ]] <- M + A %*% X %*% B
  }
  Y

}
