mn.reg <- function(Y, X, tol = 1e-6, max_iter = 1000) {

  tic <- proc.time()
  k <- length(Y)
  dm <- dim(Y[[1]])
  n <- dm[1]; p <- dm[2]
  q <- ncol(X[[1]])

  A0 <- Reduce('+', lapply(X, crossprod))
  B0 <- Reduce('+', lapply(1:k, function(i) t(X[[i]]) %*% Y[[i]]))
  B <- solve(A0, B0)

  E <- lapply(1:k, function(i) Y[[i]] - X[[i]] %*% B)
  U <- V <- 0
  for (i in 1:k) {
    U <- U + tcrossprod(E[[i]])
    V <- V + crossprod(E[[i]])
  }
  U <- U / (k * p)
  V <- V / (k * n)

  invU <- solve(U); invV <- solve(V)
  log_det_U <- determinant(U, log=TRUE)$modulus
  log_det_V <- determinant(V, log=TRUE)$modulus
  S <- Reduce('+', lapply(E, function(e) t(e) %*% invU %*% e))
  tr <- sum(diag(invV %*% S))
  logLik_old <- -0.5*k*n*p*log(2*pi) - 0.5*k*p*log_det_U - 0.5*k*n*log_det_V - 0.5*tr

  iter <- 1
  while (iter <= max_iter) {
    invU <- solve(U); invV <- solve(V)

    # Update B
    A <- Reduce('+', lapply(1:k, function(i) t(X[[i]]) %*% invU %*% X[[i]]))
    Bvec <- Reduce('+', lapply(1:k, function(i) t(X[[i]]) %*% invU %*% Y[[i]]))
    B <- solve(A, Bvec)

    E <- lapply(1:k, function(i) Y[[i]] - X[[i]] %*% B)

    # Update U and V
    U <- V <- 0
    for (i in 1:k) {
      U <- U + E[[i]] %*% invV %*% t(E[[i]])
      V <- V + t(E[[i]]) %*% invU %*% E[[i]]
    }
    U <- U / (k * p)
    V <- V / (k * n)

    # New log-likelihood
    invU <- solve(U); invV <- solve(V)
    log_det_U <- determinant(U, log=TRUE)$modulus
    log_det_V <- determinant(V, log=TRUE)$modulus
    S <- Reduce('+', lapply(E, function(e) t(e) %*% invU %*% e))
    tr <- sum(diag(invV %*% S))
    logLik_new <- -0.5*k*n*p*log(2*pi) - 0.5*k*p*log_det_U - 0.5*k*n*log_det_V - 0.5*tr

    if (abs(logLik_new - logLik_old) < tol)  break
    logLik_old <- logLik_new
    iter <- iter + 1
  }

  # Add names to B
  if (is.null(colnames(X[[1]]))) rownames(B) <- paste0("X", 1:q)
  else rownames(B) <- colnames(X[[1]])
  if (is.null(colnames(Y[[1]]))) colnames(B) <- paste0("Y", 1:p)
  else colnames(B) <- colnames(Y[[1]])

  # Fitted values (list of matrices)
  fitted <- lapply(1:k, function(i) X[[i]] %*% B)

  runtime <- proc.time() - tic
  out <- list(runtime = runtime, iters = iter, B = B, U = U, V = V,
              fitted = fitted, logLik = logLik_new)
  class(out) <- "mn.reg_rowwise"
  return(out)
}

# Predict method for rowwise model
predict.mn.reg_rowwise <- function(object, newX, ...) {
  # newX : list of matrices (each n x q) for new observations
  if (!inherits(object, "mn.reg_rowwise")) stop("Not a rowwise mn.reg object")
  lapply(newX, function(x) x %*% object$B)
}
