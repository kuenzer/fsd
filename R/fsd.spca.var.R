#' Calculate the Variance explained by Spectral Principal Components
#'
#' This function estimates the (theoretical) fraction of the variance explained
#' by each spectral principal component for some spatial functional data.
#' @param F the spectral density.
#' @return the fraction of the variance explained by each PC.
#' @seealso \code{\link{fsd.spca}}
#' @details To ensure accuracy of numerical integration, the frequencies of
#'   \code{F} should be a suitably dense grid. In order for the estimations to
#'   conform with the actual variance explained (1 - NMSE), the tuning
#'   parameters in the estimation of the spectral density need to be chosen
#'   carefully.
#' @keywords fsd
#' @export
#' @examples
#' \dontrun{
#' fsd.spca.var(F)
#' }

fsd.spca.var = function (F)
{
  if (!inherits(F, "fsd.freqdom"))
    stop("F needs to be a spectral density operator")

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

  # Integrate
  E.int = apply(pmax(Re(F.E), 0), 1, sum)
  names(E.int) = paste("PC", 1:dim(F.E)[1])

  E.int/sum(E.int)
}
