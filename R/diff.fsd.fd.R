#' Calculate the differences for a Functional Spatial Data Object
#'
#' This function is used to calculate the differences of functional data on a
#' spatial grid.
#' @param x the functional spatial data.
#' @param dim the dimension with respect to which the differences are to be
#'   computed.
#' @param lag the lag.
#' @param differences the order of the differences.
#' @param ... currently unused
#' @seealso \code{\link{fsd.fd}}
#' @keywords fsd
#' @export
#' @examples
#' data("temp")
#' plot(diff(temp, dim = 1))

diff.fsd.fd = function (x, lag = 1, differences = 1, ..., dim = 1)
{
  if (!inherits(x, "fd"))
    stop("x needs to be functional data")

  r = length(dim(x$coefs)) - 1

  if (dim > r || dim < 1)
    stop("Please select a valid dimension")

  if (dim(x$coefs)[1+dim] <= 1 + differences * lag)
    stop("Not enough data for the computation of the differences")

  pind = c(1+1:dim, 1)
  if (dim < r)
    pind = c(pind,  1+(dim+1):r)

  fdnames = x$fdnames
  if (length(fdnames[[1 + dim]]) > lag)
    fdnames[[1 + dim]] = fdnames[[1 + dim]][-(1:lag)]
  names(fdnames)[1 + dim] = paste0(names(fdnames)[1 + dim],
                                   " (differences of order ", differences,
                                   " with lag ", lag, ")")

  diffcoef = apply(x$coefs, -(1 + dim), diff,
                   lag = lag, differences = differences)

  fsd.fd(coef = aperm(diffcoef, pind),
         basisobj = x$basis, fdnames = fdnames)
}
