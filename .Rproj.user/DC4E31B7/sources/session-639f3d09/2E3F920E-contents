mn.equaltest <- function(Y1, Y2, tol = 1e-6) {

  Y <- c(Y1, Y2)
  m0 <- MN::mn.mle(Y, tol = tol)
  fit1 <- MN::mn.mle(Y1, tol = tol)
  fit2 <- MN::mn.mle(Y2, tol = tol)
  statistic <-  2 * (fit1$loglik + fit2$loglik - m0$loglik )
  n <- dim( Y1[[ 1 ]] )[1]   ;   p <- dim( Y1[[ 1 ]] )[2]
  dof <- n * p + n * (n + 1)/2 + p * (p + 1)/2

  p.value <- pchisq(statistic, dof, lower.tail = FALSE)
  names(statistic) <- "LR-test statistic"
  parameter <- dof     ;   names(parameter) <- "df"
  alternative <- "The parameters are not equal"
  method <- "Log-likelihood ratio test for equality of parameters"
  data.name <- c("data ", "groups")
  result <- list( statistic = statistic, parameter = parameter, p.value = p.value,
                  alternative = alternative, method = method, data.name = data.name )
  class(result) <- "htest"
  return(result)

}
