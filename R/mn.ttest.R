mn.ttest <- function(Y1, Y2) {
  Y <- c(Y1, Y2)
  X1 <- X2 <- list
  dm <- dim(Y[[ 1 ]])
  for ( i in 1:length(Y1) )  X1[[ i ]] <- matrix(0, nrow = dm[1], ncol = 1 )
  for ( i in 1:length(Y2) )  X2[[ i ]] <- matrix(1, nrow = dm[1], ncol = 1 )
  X <- c(X1, X2)
  stat <- 2 * ( MN::mn.reg(Y, X)$loglik - MN::mn.reg_null(Y)$loglik )
  pvalue <- pchisq(stat, dm[2], lower.tail = FALSE)
}
