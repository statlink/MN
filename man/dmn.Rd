\name{Density of the matrix normal distribution}
\alias{dmn}
\title{
Density of the matrix normal distribution
}

\description{
Density of the matrix normal distribution.
}

\usage{
dmn(X, M, U, V, logged = FALSE)
}

\arguments{
\item{X}{
A list with k elements, k matrices of dimension \eqn{n \ times p} each. In the case of one matrix only, this may be given as a numerical matrix and not as an element in a list.
}
\item{M}{
The mean matrix of the distribution, a numerical matrix of dimensions \eqn{n \times p}.
}
\item{U}{
The covariance matrix associated with the rows, a numerical matrix of dimensions \eqn{n \times n}.
}
\item{V}{
The covariance matrix associated with the columns, a numerical matrix of dimensions \eqn{p \times p}.
}
\item{logged}{
Should the logarithm of the density be computed?
}
}

%\details{
%The constraint is that all beta coefficients (including the constant) are positive.
%}

\value{
A numeric vector with the (logged) density values.
}

\references{
https://en.wikipedia.org/wiki/Matrix_normal_distribution#Definition

Pocuca, N., Gallaugher, M. P., Clark, K. M., & McNicholas, P. D. (2019). Assessing and Visualizing Matrix Variate Normality. arXiv:1910.02859.
}

\author{
Omar Alzeley.

R implementation and documentation: Omar Alzeley \email{oazeley@uqu.edu.sa}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{rmn}, \link{mn.mle}, \link{ddplot}
}
}

\examples{
M <- as.matrix(iris[1:8, 1:4])
U <- cov( matrix( rnorm(100 * 8), ncol = 8 ) )
V <- cov( iris[1:50, 1:4] )
X <- rmn(10, M, U, V)
dmn(X, M, U, V, TRUE)
}
