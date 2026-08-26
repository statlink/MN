\name{Distance-Distance Plot}
\alias{ddplot}
\title{
Distance-Distance Plot
}

\description{
Distance-Distance Plot
}

\usage{
ddplot(X, M, U, V)
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
}

\details{
The distance-distance plot is produced. This is a scatter plot of the Mahalanobis distances computed using the estimated parameters from the multivariate normal and matrix normal distribution. See Pocuca et al. (2019) for more details.
}

\value{
A scatter plot of the Mahalanobis distances.
}

\references{
Pocuca N., Gallaugher M. P., Clark K. M. & McNicholas P. D. (2019). Assessing and Visualizing Matrix Variate Normality. arXiv:1910.02859.
}

\author{
Michail Tsagris and Omar Alzeley.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr} and Omar Alzeley \email{oazeley@uqu.edu.sa}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{rmn}, \link{mn.mle}, \link{dmn}, \link{ddkstest}
}
}

\examples{
M <- as.matrix(iris[1:8, 1:4])
U <- cov( matrix( rnorm(100 * 8), ncol = 8 ) )
V <- cov( iris[1:50, 1:4] )
X <- rmn(100, M, U, V)
ddplot(X, M, U, V)
}
