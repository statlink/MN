\name{Random matrices simulation from the matrix normal distribution}
\alias{rmn}
\title{
Random matrices simulation from the matrix normal distribution
}

\description{
Random matrices simulation from the matrix normal distribution.
}

\usage{
rmn(k, M, U, V)
}

\arguments{
\item{k}{
The sample size, the number of matrices to simulate.
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
}

%\details{
%This is described in Chapter 8.2 of Hansen (2019). The idea is to inimise the sum of squares of the residuals under the constraint \eqn{R^T \beta = c}. As mentioned above, be careful with the input you give in the x matrix and the R vector.
%}

\value{
A list with k elements, k matrices of dimension \eqn{n \ times p} each. These are the random matrices drawn from a matrix normal distribution.
}

\references{
https://en.wikipedia.org/wiki/Matrix_normal_distribution#Definition
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{dmn}, \link{mn.mle}, \link{ddplot}
}
}

\examples{
M <- as.matrix(iris[1:8, 1:4])
U <- cov( matrix( rnorm(100 * 8), ncol = 8 ) )
V <- cov( iris[1:50, 1:4] )
X <- rmn(10, M, U, V)
}
