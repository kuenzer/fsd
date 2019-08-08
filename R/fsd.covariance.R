#' Estimate autocovariance operators
#'
#' This function is used to estimate a covariance operator of stationary
#' functional spatial data with respect to some lag h.
#' @param Y the functional spatial data. Either an fd object or an array.
#' @param h the spatial lag with respect to which the autocovariance is to be
#'   estimated.
#' @param centered a boolean indicating whether or not the data is already
#'   centered.
#' @param na.rm whether or not missing values should be ignored.
#' @param unbiased whether or not to scale the estimator to be unbiased. FALSE
#'   by default. Is not used if na.rm is TRUE.
#' @return the covariance matrix.
#' @keywords fsd
#' @export
#' @examples
#' \dontrun{
#' fsd.covariance(Y, h)
#' }

fsd.covariance = function (Y = NULL, h, centered = FALSE, na.rm = FALSE,
                           unbiased = FALSE)
{
  if (inherits(Y, "fd"))
    Y = Y$coefs

  r = length(dim(Y)) - 1
  n = dim(Y)[1+1:r]
  d = dim(Y)[1]

  if (any(n - abs(h) <= 0))
    return(NULL)

  if (!centered)
    Y = Y - apply(Y, 1, mean)

  if (na.rm) {
    Yna = apply(is.na(Y), 1+1:r, any)
    Y[1:dim(Y)[1] + dim(Y)[1] * rep(which(Yna) - 1, each = dim(Y)[1])] = 0
  }

  ind2 = mapply(function(ni, hi){max(-hi, 0) + 1:(ni - abs(hi))}, n, h,
                SIMPLIFY = FALSE)
  ind1 = mapply("+", ind2, h, SIMPLIFY = FALSE)

  if (any(sapply(ind1, length) > 1)) {
    Ch = crossprod(do.call(what = apply, args = c(X = list(Y), MARGIN = 1,
                                          FUN = '[', ind1)),
                   do.call(what = apply, args = c(X = list(Y), MARGIN = 1,
                                          FUN = '[', ind2)))
  } else { # if the mean is taken over 1 pair, make sure multiplication works
    Ch = tcrossprod(do.call(what = apply, args = c(X = list(Y), MARGIN = 1,
                                        FUN = '[', ind1)),
                    do.call(what = apply, args = c(X = list(Y), MARGIN = 1,
                                        FUN = '[', ind2)))
  }

  if (na.rm)
    Ch / sum( !(do.call('[', c(list(Yna), ind1)) |
                  do.call('[', c(list(Yna), ind2))) )
  else if (unbiased)
    Ch / prod(n - abs(h))
  else
    Ch / prod(n)
}
