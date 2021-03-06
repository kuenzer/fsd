% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fsd.spca.cov.R
\name{fsd.spca.cov}
\alias{fsd.spca.cov}
\title{Calculate the Covariance of the Spectral Principal Component Scores}
\usage{
fsd.spca.cov(F, L = 3)
}
\arguments{
\item{F}{the spectral density.}

\item{L}{the maximum lag to compute the covariance of the scores}
}
\value{
A list with components \item{laglist}{the list of lags.} \item{cov}{a
  matrix with the autocovariances of the principal component scores.}
}
\description{
This function estimates the (theoretical) autocovariance for each spectral
principal component score for some spatial functional data.
}
\details{
To ensure accuracy of numerical integration during the Fourier
  transform, the frequencies of \code{F} should be a suitably dense grid.
}
\examples{
\dontrun{
fsd.spca.cov(F, L)
}
}
\seealso{
\code{\link{fsd.spca}}, \code{\link{fsd.spca.var}}
}
\keyword{fsd}
