#' Compute the centered version of spatial functional data
#'
#' This function is used to center spatial functional data.
#' @param X the functional spatial data.
#' @param margins the margins for which the centering should be done separately
#'   in every entry, in order to remove a trend.
#' @param na.rm whether or not missing values should be ignored.
#' @return the centered data.
#' @seealso \code{\link{mean.fsd.fd}}
#' @keywords fsd
#' @export
#' @examples
#' data("temp")
#'
#' # blindly center the data
#' T.centered = center.fsd.fd(temp)
#' # there is a clear structure in the empirical means
#' plot(mean(T.centered, margins = 1:2))
#' 
#' # take the latitude and the longitude as baseline
#' T.centered2 = center.fsd.fd(temp, margins = 1:2)
#' # random noise in the empirical means
#' plot(mean(T.centered2, margins = 1:2))
#' # this is the baseline
#' plot(mean(temp, margins = 1:2))

center.fsd.fd = function (X, margins = NULL, na.rm = FALSE)
{
  if (!inherits(X, "fd"))
    stop("X needs to be functional data")
  if (is.null(margins))
    return(fsd.fd(X$coefs - as.vector(mean.fsd.fd(X, na.rm = na.rm)$coefs),
                  basisobj = X$basis, fdnames = X$fdnames))
  else {
    if (any(margins > length(dim(X$coefs))-1) || any(margins < 1))
      stop("Please provide valid margins")

    r = length(dim(X$coefs)) - 1
    permut = c(1, margins +1, (1:r)[-margins] +1)
    Y = aperm(X$coefs, permut)
    Y.mean = apply(Y, seq(1+length(margins)), mean, na.rm = na.rm)
    Y = aperm(Y-as.vector(Y.mean), order(permut))

    return(fsd.fd(coef = Y,
                  basisobj = X$basis, fdnames = X$fdnames))
  }
}
