dmn <- function(X, M, U, V, logged = FALSE) {

  invU <- solve(U)  ;  invV <- solve(V)

  if ( !is.list(X) ) {
    X <- as.matrix(X)
    dm <- dim(X)
    n <- dm[1]  ;  p <- dm[2]
    Z <- X - M
    a <- invV %*% t(Z) %*% invU %*% Z
    f <-  -0.5 * sum( diag(a) )

  } else {

    k <- length(X)
    dm <- dim(X[[ 1 ]])
    n <- dm[1]  ;  p <- dm[2]
    f <- numeric(k)

    Z <- lapply(X, M, FUN = "-")
    for ( i in 1:k ) {
      a <- invV %*% t( Z[[ i ]] ) %*% invU %*% Z[[ i ]]
      f[i] <-  -0.5 * sum( diag(a) )
    }
  }

  con <- 0.5 * n * p * log(2 * pi) + 0.5 * n * log( det(V) ) +
    0.5 * p * log( det(U) )
  f <- f - con
  if ( !logged )  f <- exp(f)
  f
}

