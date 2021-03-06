% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fsd.plot.data.R
\name{fsd.plot.data}
\alias{fsd.plot.data}
\title{Plot functional spatial data}
\usage{
fsd.plot.data(
  X,
  basisobj = NULL,
  gridsize = NULL,
  set.dims = NULL,
  xlim = NULL,
  ylim = NULL,
  main = NULL,
  xlab = NULL,
  ylab = NULL
)
}
\arguments{
\item{X}{the functional spatial data. Either an fd object or an array, or a
list thereof.}

\item{basisobj}{the basis of the functional data.}

\item{gridsize}{the number of grid points on which to evaluate the functional
data.}

\item{set.dims}{the coordinates to be fixed. A vector of dimension r. NaN
represents not fixed components.}

\item{xlim}{limits for the coordinates on the x-axis.}

\item{ylim}{limits for the coordinates on the y-axis.}

\item{main}{title of the plot.}

\item{xlab}{label for the x-axis.}

\item{ylab}{label for the y-axis.}
}
\description{
This function is used to plot functional spatial data on a grid of dimension
r,  along a two-dimensional plane of indices.
}
\examples{
data("temp")
# plot the data for the third year in the dataset (i.e. 1981)
fsd.plot.data(temp, set.dims = c(NaN, NaN, 3))
}
\seealso{
\code{\link{plot.fsd.fd}}, \code{\link{fsd.plot.filters}},
  \code{\link{fsd.plot.covariance}}
}
\keyword{fsd}
