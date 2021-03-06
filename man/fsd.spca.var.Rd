% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fsd.spca.var.R
\name{fsd.spca.var}
\alias{fsd.spca.var}
\title{Calculate the Variance explained by Spectral Principal Components}
\usage{
fsd.spca.var(F)
}
\arguments{
\item{F}{the spectral density.}
}
\value{
the fraction of the variance explained by each PC.
}
\description{
This function estimates the (theoretical) fraction of the variance explained
by each spectral principal component for some spatial functional data.
}
\details{
To ensure accuracy of numerical integration, the frequencies of
  \code{F} should be a suitably dense grid. In order for the estimations to
  conform with the actual variance explained (1 - NMSE), the tuning
  parameters in the estimation of the spectral density need to be chosen
  carefully.
}
\examples{
\dontrun{
fsd.spca.var(F)
}
}
\seealso{
\code{\link{fsd.spca}}
}
\keyword{fsd}
