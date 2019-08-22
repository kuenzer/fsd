#' Estimate the spectral density operator
#'
#' This function is used to estimate the spectral density operator of stationary
#' functional spatial data on a spatial frequency grid.
#' @param X the functional spatial data. Either an fd object or an array.
#' @param freq the spatial frequencies with respect to which the spectral
#'   density is to be estimated. A vector or a list of vectors defining the
#'   grid.
#' @param q the size of the window of the kernel used for estimation. An integer
#'   or a vector. For any dimension, the value \code{q = 1} is equivalent to
#'   \code{q = 0} and means that no lags in this direction are taken into
#'   account.
#' @param na.rm whether or not missing values should be ignored.
#' @return the spectral density, i.e. a list consisting of \item{operators}{a
#'   list of the the spectral density operator at different spatial
#'   frequencies.} \item{freq}{a list of the spatial frequencies.}
#'   \item{basis}{the basis object of the functional data.}
#'   \item{estimation.q}{the tuning parameter \code{q} used in the estimation
#'   process.}
#' @details The spectral density is estimated employing a  Bartlett kernel on
#'   the Euclidean norm of the lags. Therefore, the lag \code{q} is the smallest
#'   lag not to be taken into consideration.
#' @keywords fsd
#' @export
#' @examples
#' \dontrun{
#' fsd.spectral.density(X, freq)
#' }

fsd.spectral.density = function (X, freq, q = NULL, na.rm = FALSE)
{
  if (inherits(X, "fd")) {
    basisobj = X$basis
    X = X$coefs
  } else
    basisobj = NULL

  d = dim(X)[1]
  r = length(dim(X)) - 1
  n = dim(X)[1+1:r]

  X = X - apply(X, 1, mean, na.rm = na.rm)

  if (is.null(q))
    q = ceiling(n^0.4)
  else if (length(q) == 1)
    q = rep(q, r)

  q = ceiling(q)

  if (length(q) != r)
    stop("Dimensions don't fit")
  if (any(q < 1))
    stop("Choose a valid q")

  hlist = unfold(lapply(q - 1, function(w) {-w:w}))

  w = sapply(hlist, function(u) {max(0, 1 - sqrt(sum((u/q)^2)))})

  hlist = hlist[w > 0]
  w = w[w > 0]

  Clist = sapply(hlist, FUN = fsd.covariance, Y = X,
                 centered = TRUE, na.rm = na.rm, simplify = "array")

  Clist = Clist * rep(w, each = prod(dim(Clist)[1:2]))

  if (!is.list(freq))
    thlist = unfold(freq, r)
  else
    thlist = unfold(freq)

  hlist = lapply(hlist, "-")

  Fth = fsd.fourier(Clist, hlist, thlist)$operators

  Fobj = list(operators = Fth, freq = thlist,
              basis = basisobj, estimation.q = q)
  class(Fobj) = "fsd.freqdom"
  Fobj
}
