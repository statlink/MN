ddplot <- function(X, M, U, V) {

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

  f1 <- knots( ecdf(dm) )  ;  f2 <- knots( ecdf(d) )

  lab1 <- round( seq(min(f1), max(f1), by = 5) )
  lab2 <- round( seq(min(f1), max(f1), by = 5) )
  plot( f1, f2, xlab = "Matrix variate", ylab = "Multivariate",
        cex.lab = 1.3, cex.axis = 1.3, col = "blue", pch = 20, xaxt = "n", yaxt = "n")
  abline( v = lab1 )
  abline( h = lab2 )
  abline(a = 0, b = 1, col = "red")
  points( f1, f2, col = "blue", pch = 20)
  axis(1, at = lab1, lab1)
  axis(2, at = lab2, lab2)

}

