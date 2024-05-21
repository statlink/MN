mn.mle <- function(X) {

  tic <- proc.time()

  k <- length(X)
  dm <- dim(X[[ 1 ]])
  n <- dm[1]  ;  p <- dm[2]
  ## mean matrix
  M <- Reduce('+', X) / k
  Z <- lapply(X, M, FUN = "-")
  tZ <- lapply(Z, FUN = "t")
  ## step 1
  U <- V <- 0
  for ( i in 1:k ) {
    U <- U + tcrossprod( Z[[ i ]] )
    V <- V + crossprod( Z[[ i ]] )
  }
  U1 <- U / (k * p)
  V1 <- V / (k * n)
  ## step 2
  U2 <- V2 <- 0
  invV1 <- solve(V1)
  for ( i in 1:k )  U2 <- U2 + Z[[ i ]] %*% invV1 %*% tZ[[ i ]]
  U2 <- U2 / (k * p)
  invU2 <- solve(U2)
  for ( i in 1:k )  V2 <- V2 + tZ[[ i ]] %*% invU2 %*% Z[[ i ]]
  V2 <- V2 / (k * n)

  i <- 2
  while ( Rfast::Norm(U1 - U2) > 1e-6 | Rfast::Norm(V1 - V2) > 1e-6 ) {
    i <- i + 1
    U1 <- U2
    V1 <- V2
    U2 <- V2 <- 0
    invV1 <- solve(V1)
    for ( i in 1:k )  U2 <- U2 + Z[[ i ]] %*% invV1 %*% tZ[[ i ]]
    U2 <- U2 / (k * p)
    invU2 <- solve(U2)
    for ( i in 1:k )  V2 <- V2 + tZ[[ i ]] %*% invU2 %*% Z[[ i ]]
    V2 <- V2 / (k * n)
  }

  runtime <- proc.time() - tic

  list(runtime = runtime, iters = i, M = M, U = U2, V = V2)
}





