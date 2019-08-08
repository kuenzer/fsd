#' Calculate the Covariance of the Spectral Principal Component Scores
#'
#' This function estimates the (theoretical) autocovariance for each spectral
#' principal component score for some spatial functional data.
#' @param F the spectral density.
#' @param L the maximum lag to compute the covariance of the scores
#' @return A list with components \item{laglist}{the list of lags.} \item{cov}{a
#'   matrix with the autocovariances of the principal component scores.}
#' @seealso \code{\link{fsd.spca}}, \code{\link{fsd.spca.var}}
#' @details To ensure accuracy of numerical integration during the Fourier
#'   transform, the frequencies of \code{F} should be a suitably dense grid.
#' @keywords fsd
#' @export
#' @examples
#' \dontrun{
#' fsd.spca.cov(F, L)
#' }

fsd.spca.cov = function (F, L = 3)
{
  if (!inherits(F, "fsd.freqdom"))
    stop("F needs to be a spectral density operator")

  if (length(L) == 1)
    L = rep(L, length(F$freq[[1]]))

  if (length(L) != length(F$freq[[1]]))
    stop("Dimensions incompatible")

  # Change of base
  if (!is.null(F$basis)) {
    B = inprod(F$basis, F$basis)
  } else {
    B = diag(dim(F$operators[[1]])[1])
  }
  B.root = eigen(B)$vectors %*% diag(sqrt(eigen(B)$values)) %*%
    t(eigen(B)$vectors)
  for (i in 1:length(F$operators)) {
    F$operators[[i]] = B.root %*% F$operators[[i]] %*% B.root
  }

  # Compute eigenvalues
  F.E = sapply(F$operators,
               function(w) {eigen(w, symmetric = TRUE,
                                  only.values = TRUE)$values},
               simplify="array")

  F.E = array(F.E, dim = c(1, dim(F.E)))

  laglist = unfold(lapply(L, function(h) {-h:h}))

  filters = fsd.fourier.inverse(F.E, F$freq, laglist)

  covs = array(unlist(filters$operators),
               dim = c(dim(F.E)[2], length(filters$operators)),
               dimnames = list(paste("PC", 1:dim(F.E)[2]), NULL))

  list(laglist = laglist, cov = covs)
}
