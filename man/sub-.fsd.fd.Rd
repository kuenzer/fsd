% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fsd.fd.R
\name{[.fsd.fd}
\alias{[.fsd.fd}
\title{Subsetting a functional spatial data object}
\usage{
\method{[}{fsd.fd}(fsdobj, ..., drop = TRUE)
}
\arguments{
\item{fsdobj}{a functional spatial data object}

\item{...}{the indices to take}

\item{drop}{whether to flatten the spatial grid by dropping dimensions of
length one.}
}
\description{
Subsetting works as with conventional arrays. By default, dimensions can be 
dropped.
}
\examples{
data("temp")
plot(temp[1:2, 1:2, 1, drop = TRUE])
}
\keyword{fsd}
