#' Perform Spectral Principal Components Analysis on Spatial Functional Data
#'
#' This function performs spectral PCA on functional spatial data on a grid of
#' dimension r.
#' @param X the functional spatial data. Either an fd object or an array.
#' @param freq.res the resolution for the computation of the spectral density.
#' @param Npc the number of principal components to be computed.
#' @param L the maximum lag for the filters. An integer or vector of integers.
#' @param q a tuning parameter for the estimation of the the spectral density
#'   operator. An integer or vector of integers.
#' @param na.ignore whether to ignore missing data points in the computation of
#'   the scores.
#' @param return.F a boolean indicating whether to return the spectral density.
#' @param only.filters a boolean indicating whether to compute only the filters,
#'   leaving out the scores.
#' @return A list with components \item{F}{the spectral density.}
#'   \item{tuning.params}{a list of the used tuning parameters \code{freq.res},
#'   \code{q} and \code{Lmax}.} \item{filters}{the SPC filters.}
#'   \item{scores}{the SPC scores.} \item{var}{the theoretical fractions of
#'   variance explained by each PC.} \item{X.mean}{the mean of X.}
#' @seealso \code{\link{fsd.spca.inverse}}, \code{\link{fsd.spectral.density}}
#' @details This function can be used to compute spectral PCA. By setting
#'   \code{q = 0}, it can also be used to perform static PCA.
#'
#'   Setting a higher value for \code{freq.res} increases the runtime
#'   significantly, but yields a more accurate result because of reduced
#'   integration errors.
#' @keywords fsd
#' @export
#' @examples
#' \dontrun{
#' fsd.spca(X)
#' }

fsd.spca = function (X, freq.res = 100, Npc = 3, L = 3, q = NULL,
                     na.ignore = TRUE, return.F = FALSE, only.filters = FALSE)
{
  # Check parameters
  if (any(L < 0))
    stop("Please choose a valid parameter L")
  if (Npc <= 0)
    stop("Please choose a valid parameter Npc")
  if (any(freq.res < 0))
    warning("Please choose a valid frequency resolution")

  # Check data type and center the data
  if (inherits(X, "fd")) {
    basisobj = X$basis
    meanX = mean.fsd.fd(X, na.rm = TRUE)
    X = X$coefs
    X = X - as.vector(meanX$coefs)
  } else {
    basisobj = NULL
    meanX = apply(X, 1, mean, na.rm = TRUE)
    X = X - meanX
  }

  na.rm = FALSE
  if (anyNA(X)) {
    na.rm = TRUE
    if (na.ignore)
      warning("The data contains missing points which will be ignored ",
              "in the computation of the SPC filters and scores.")
    else
      warning("The data contains missing points which will be ignored ",
              "in the computation of the SPC filters but not the scores.")
  }
  na.ignore = na.ignore & na.rm

  # Adjust the tuning parameters
  L = floor(L)
  if (length(L) == 1)
    L = rep(L, length(dim(X)) - 1)

  if (length(freq.res) == 1)
    freq.res = rep(freq.res, length(dim(X)) - 1)

  if (!is.null(q)) {
    q = pmax(q, 1)
    q = floor(q)
    freq.res = freq.res * (q > 1)
    L = L * (q > 1)
  }

  if (any(freq.res < 3 * L))
    warning("The chosen frequency resolution is very low.\n",
            "Significant errors during the integration may occur.")

  # Compute the spectral density and the PC filters
  freq = lapply(freq.res,
                function(fr) {if (fr > 0) -(fr - 1):fr / fr * pi else 0})

  dims = dim(X)[1]
  Npc = min(Npc, dims)

  F = fsd.spectral.density(X, freq, q = q, na.rm = na.rm)
  F$basis = basisobj

  filters = fsd.spca.filters(F, Npc, L)

  # Compute the PC scores
  scores = NULL
  if (!only.filters)
    scores = fsd.spca.scores(X, filters, na.ignore = na.ignore)

  Fvar = filters$var
  filters$var = NULL

  qF = F$estimation.q

  if (!return.F)
    F = NULL

  na.status = ifelse(na.rm, ifelse(na.ignore, "ignored", "removed"), "none")

  list(F = F, tuning.params = list(freq.res = freq.res,
                                   estimation.q = qF, Lmax = L, na = na.status),
       filters = filters, scores = scores, var = Fvar, X.mean = meanX)
}
