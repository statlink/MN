mn.ttest <- function(Y1, Y2, tol = 1e-6, maxiter = 100) {
  Y <- c(Y1, Y2)
  X1 <- Y1   ;   X2 <- Y2
  dm <- dim( Y1[[ 1 ]] )
  for ( i in 1:length(Y1) )  X1[[ i ]] <- matrix(0, nrow = dm[1], ncol = 1 )
  for ( i in 1:length(Y2) )  X2[[ i ]] <- matrix(1, nrow = dm[1], ncol = 1 )
  X <- c(X1, X2)
  statistic <- 2 * ( MN::mn.reg(Y, X, tol = tol, maxiter = maxiter)$loglik -
                     MN::mn.reg_null(Y, tol = tol, maxiter = maxiter)$loglik )
  p.value <- pchisq(statistic, dm[2], lower.tail = FALSE)
  names(statistic) <- "LR-test statistic"
  parameter <- dm[2]     ;   names(parameter) <- "df"
  alternative <- "The mean matrices differ"
  method <- "Log-likelihood ratio test for equality of mean matrices"
  data.name <- c("data ", "groups")
  result <- list( statistic = statistic, parameter = parameter, p.value = p.value,
                  alternative = alternative, method = method, data.name = data.name )
  class(result) <- "htest"
  return(result)
}
