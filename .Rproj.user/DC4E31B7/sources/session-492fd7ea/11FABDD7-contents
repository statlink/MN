mn.reg <- function(Y, X, tol = 1e-6, max_iter = 100, Xnew = NULL) {

  tic <- proc.time()
  k <- length(Y)
  dm <- dim( Y[[ 1 ]] )
  n <- dm[1]   ;   p <- dm[2]
  dx <- ncol( X[[ 1 ]] )[2]

  A0 <- Reduce( '+', lapply(X, crossprod) )
  B0 <- Reduce( '+', lapply( 1:k, function(i)  t( X[[ i ]]) %*% Y[[ i ]] ) )
  B <- solve(A0, B0)
  E <- lapply( 1:k, function(i)  Y[[ i ]] - X[[ i ]] %*% B )
  U <- V <- 0
  for ( i in 1:k ) {
    U <- U + tcrossprod( E[[ i ]] )
    V <- V + crossprod( E[[ i ]] )
  }
  U <- U / (k * p)
  V <- V / (k * n)

  invU <- solve(U)   ;   invV <- solve(V)
  ldU <- as.numeric( determinant(U, logarithm = TRUE)$modulus )
  ldV <- as.numeric( determinant(V, logarithm = TRUE)$modulus )
  S <- Reduce( '+', lapply( E, function(e)  t(e) %*% invU %*% e ) )
  tr <- sum( diag( invV %*% S ) )
  loglik_old <- -0.5 * k * p * ldU - 0.5 * k * n * ldV - 0.5 * tr
  con <-  -0.5 * k * n * p * log(2 * pi)

  iter <- 1
  while ( iter <= max_iter ) {
    invU <- solve(U)   ;   invV <- solve(V)
    A <- Reduce( '+', lapply(1:k, function(i) t( X[[ i ]] ) %*% invU %*% X[[ i ]] ) )
    C <- Reduce( '+', lapply(1:k, function(i) t( X[[ i ]] ) %*% invU %*% Y[[ i ]] ) )
    B <- solve(A, C)
    E <- lapply( 1:k, function(i)  Y[[ i ]] - X[[ i ]] %*% B )

    U <- 0
    for ( i in 1:k )  U <- U + E[[ i ]] %*% invV %*% t(E[[ i ]])
    U <- U / (k * p)

    invU_new <- solve(U)
    V <- 0
    for ( i in 1:k )  V <- V + t( E[[ i ]] ) %*% invU_new %*% E[[ i ]]
    V <- V / (k * n)

    invU <- solve(U)   ;   invV <- solve(V)
    ldU <- as.numeric( determinant(U, logarithm = TRUE)$modulus )
    ldV <- as.numeric( determinant(V, logarithm = TRUE)$modulus )
    S <- Reduce( '+', lapply( E, function(e)  t(e) %*% invU %*% e ) )
    tr <- sum( diag(invV %*% S) )
    loglik_new <- -0.5 * k * p * ldU - 0.5 * k * n * ldV - 0.5 * tr

    # Relative change convergence
    if (abs( loglik_new - loglik_old) / abs(loglik_old) < tol )  break
    loglik_old <- loglik_new
    iter <- iter + 1
  }

  if ( is.null( colnames( X[[ 1 ]] ) ) )  rownames(B) <- paste0("X", 1:dx)
  else  rownames(B) <- colnames( X[[ 1 ]] )
  if ( is.null( colnames( Y[[ 1 ]] ) ) )  colnames(B) <- paste0("Y", 1:p)
  else  colnames(B) <- colnames( Y[[ 1 ]] )

  est <- NULL
  if ( !is.null(Xnew) )  est <- lapply( Xnew, function(x) x %*% B )

  runtime <- proc.time() - tic

  list( runtime = runtime, iters = iter, loglik = con + loglik_new, B = B,
    U = U, V = V, est = est )
}
