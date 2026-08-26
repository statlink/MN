\name{Log-likelihood ratio test for equality of mean matrices}
\alias{mn.ttest}
\title{
Log-likelihood ratio test for equality of mean matrices
}

\description{
Log-likelihood ratio test for equality of mean matrices.
}

\usage{
mn.ttest(Y1, Y2, tol = 1e-6, maxiter = 100)
}

\arguments{
\item{Y1}{
The observations of the first sample. A list with k1 elements (k1 is the sample size),
k1 matrices of dimension \eqn{n \times p} each.
}
\item{Y2}{
The observations of the first sample. A list with k2 elements (k2 is the sample size),
k2 matrices of dimension \eqn{n \times p} each.
}
\item{tol}{
The tolerance value to terminate the iterative algorithm.
}
\item{maxiter}{
The maximum number of iterations allowed.
}
}

\details{
A log-likelihood ratio test comparing the equality of mean matrices, \eqn{M_1=M_2} is performed,
and the alternative is that the mean matrices are not equal.
}

\value{
This is an "htest"class object. Thus it returns a list including:
\item{statistic}{
The test statistic value.
}
\item{parameter}{
The degrees of freedom of the test.
}
\item{p.value}{
The p-value of the test.
}
\item{alternative}{
A character with the alternative hypothesis.
}
\item{method}{
A character with the test used.
}
\item{data.name}{
A character vector with two elements.
}
}

\references{
https://en.wikipedia.org/wiki/Matrix_normal_distribution#Definition

Pocuca N., Gallaugher M. P., Clark K. M. & McNicholas P. D. (2019).
Assessing and Visualizing Matrix Variate Normality. arXiv:1910.02859.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{mn.mle}, \link{mn.reg}
}
}

\examples{
M <- as.matrix(iris[1:8, 1:4])
U <- cov( matrix( rnorm(100 * 8), ncol = 8 ) )
V <- cov( iris[1:50, 1:4] )
Y1 <- MN::rmn(150, M, U, V)
Y2 <- MN::rmn(100, M, U, V)
mn.ttest(Y1, Y2)
}
