% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{temp}
\alias{temp}
\title{Temperature Data in Wyoming}
\format{An object of class \code{"fsd.fd"}.}
\source{
\url{https://www.esrl.noaa.gov/psd/data/gridded/data.cpc.globaltemp.html}
}
\usage{
data(temp)
}
\description{
CPC daily Temperature data of mean temperatures in Wyoming
projected onto a B-spline basis of 36 functions for each year,
from 1979 to 2017 and with a spatial resolution of 0.5 degrees.
}
\details{
The gridded data has been obtained by Shepard's method.
}
\examples{
data(temp)
plot(temp)
}
\references{
ESRL Physical Sciences Division
(\url{https://www.esrl.noaa.gov/psd/})
}
\keyword{datasets}
