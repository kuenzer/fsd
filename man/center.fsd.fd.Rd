% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/center.fsd.fd.R
\name{center.fsd.fd}
\alias{center.fsd.fd}
\title{Compute the centered version of spatial functional data}
\usage{
center.fsd.fd(X, margins = NULL, na.rm = FALSE)
}
\arguments{
\item{X}{the functional spatial data.}

\item{margins}{the margins for which the centering should be done separately
in every entry, in order to remove a trend.}

\item{na.rm}{whether or not missing values should be ignored.}
}
\value{
the centered data.
}
\description{
This function is used to center spatial functional data.
}
\examples{
data("temp")

# blindly center the data
T.centered = center.fsd.fd(temp)
# there is a clear structure in the empirical means
plot(mean(T.centered, margins = 1:2))

# take the latitude and the longitude as baseline
T.centered2 = center.fsd.fd(temp, margins = 1:2)
# random noise in the empirical means
plot(mean(T.centered2, margins = 1:2))
# this is the baseline
plot(mean(temp, margins = 1:2))
}
\seealso{
\code{\link{mean.fsd.fd}}
}
\keyword{fsd}
