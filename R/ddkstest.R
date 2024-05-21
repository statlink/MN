ddkstest <- function(X, M, U, V, alpha = 0.05) {

  invU <- solve(U)  ;  invV <- solve(V)
  k <- length(X)
  dm <- dim(X[[ 1 ]])
  n <- dm[1]  ;  p <- dm[2]
  dm <- numeric(k)

  x <- matrix( nrow = k, ncol = n * p )

  Z <- lapply(X, M, FUN = "-")
  for ( i in 1:k ) {
    a <- invU %*% Z[[ i ]] %*% invV %*% t( Z[[ i ]] )
    dm[i] <- sum( diag(a) )
    x[i, ] <- as.vector(X[[ i ]])
  }

  m <- as.vector(M)
  s <- crossprod( Rfast::eachrow(x, m, oper = "-") ) / (k - 1)
  d <- Rfast::mahala(x, m, s)

  f1 <- ecdf(dm)  ;  f2 <- ecdf(d)
  q1 <- f1(dm)  ;  q2 <- f2(d)
  res <- max( abs(q1 - q2) ) > sqrt( -0.5 * log(alpha) * 2 * k/k^2 )
  if (res) {
    paste("Reject")
  } else paste ("Not reject")
}

