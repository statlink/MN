\name{Maximum likelihood estimation of the the matrix normal distribution}
\alias{mn.mle}
\title{
Maximum likelihood estimation of the the matrix normal distribution
}

\description{
Maximum likelihood estimation of the the matrix normal distribution.
}

\usage{
mn.mle(X)
}

\arguments{
\item{X}{
A list with k elements (k is the sample size), k matrices of dimension \eqn{n \times p} each.
}
}

%\details{
%The constraint is that all beta coefficients (including the constant) are positive.
%}

\value{
A list including:
\item{runtime}{
The runtime required for the whole fitting procedure.
}
\item{iters}{
The number of iterations required for the estimation of the U and V matrices.
}
\item{M}{
The estimated mean matrix of the distribution, a numerical matrix of dimensions \eqn{n \times p}.
}
\item{U}{
The estimated covariance matrix associated with the rows, a numerical matrix of dimensions \eqn{n \times n}.
}
\item{V}{
The estimated covariance matrix associated with the columns, a numerical matrix of dimensions \eqn{p \times p}.
}
}

\references{
https://en.wikipedia.org/wiki/Matrix_normal_distribution#Definition

Pocuca N., Gallaugher M. P., Clark K. M. & McNicholas P. D. (2019). Assessing and Visualizing Matrix Variate Normality. arXiv:1910.02859.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{dmn}, \link{rmn}, \link{ddplot}
}
}

\examples{
M <- as.matrix(iris[1:8, 1:4])
U <- cov( matrix( rnorm(100 * 8), ncol = 8 ) )
V <- cov( iris[1:50, 1:4] )
X <- rmn(200, M, U, V)
mod <- mn.mle(X)
}
