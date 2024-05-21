\name{Kolmogorov-Smirnov test for matrix normality}
\alias{ddkstest}
\title{
Kolmogorov-Smirnov test for matrix normality
}

\description{
Kolmogorov-Smirnov test for matrix normality
}

\usage{
ddkstest(X, M, U, V, alpha = 0.05)
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
\item{alpha}{
The significance level for the test, set by default equal to 0.05.
}
}

\details{
The Kolmogorov-Smirnov test for matrix normality is performed. See Pocuca (2019) for more details.
}

\value{
A message. If the Kronecker product covariance structure is not present, the message reads "Reject" and "Not reject otherwise".
}

\references{
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
\code{\link{rmn}, \link{mn.mle}, \link{dmn}, \link{ddplot}
}
}

\examples{
M <- as.matrix(iris[1:8, 1:4])
U <- cov( matrix( rnorm(100 * 8), ncol = 8 ) )
V <- cov( iris[1:50, 1:4] )
X <- rmn(200, M, U, V)
ddkstest(X, M, U, V)
}
