#' Compute the mean of spatial functional data
#'
#' This function is used to calculate the mean on functional data on a spatial
#' grid.
#' @param x the functional spatial data.
#' @param margins the margins that should not be averaged over.
#' @param trim currently unused
#' @param na.rm whether or not missing values should be ignored.
#' @param ... currently unused
#' @return the mean of x
#' @seealso \code{\link{center.fsd.fd}}
#' @keywords fsd
#' @export
#' @examples
#' data("temp")
#'
#' # the mean temperature curve
#' plot(mean(temp))
#' # slight variations with respect to the location can be seen
#' plot(mean(temp, margins = 1:2))

mean.fsd.fd = function (x, trim = 0, na.rm = FALSE, ..., margins = NULL)
{
  if (!inherits(x, "fd"))
    stop("x needs to be functional data")
  if (is.null(margins))
    fsd.fd(coef = matrix(apply(x$coefs, 1, mean, na.rm = na.rm), ncol = 1),
           basisobj = x$basis, fdnames = list(x$fdnames[[1]], "mean"))
  else {
    if (any(margins > length(dim(x$coefs))-1) || any(margins < 1))
      stop("Please provide valid margins")
    fsd.fd(coef = apply(x$coefs, c(1, 1 + margins), mean, na.rm = na.rm),
           basisobj = x$basis, fdnames = x$fdnames[c(1, 1 + margins)])
  }
}
