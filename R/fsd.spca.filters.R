#' Calculate the Spectral Principal Components Filters Spatial Functional Data
#'
#' This function calculates the SPC filters from a given spectral density
#' operator.
#' @param F the spectral density.
#' @param Npc the number of principal components to be computed.
#' @param L the maximum lag for the filters, as an integer or vector of
#'   integers.
#' @return the SPC filters.
#' @seealso \code{\link{fsd.spca}}
#' @details The eigenvectors used for the calculation of the Fourier expansion
#'   are oriented such that the sum of their coordinates over the basis of the
#'   fd object (i.e. the entries in the \code{coefs} array) is a non-negative
#'   real number. If the basis is not orthonormal, this is done with respect to
#'   a virtual orthonormal basis in the background.
#' @keywords fsd
#' @export
#' @examples
#' \dontrun{
#' fsd.spca.filters(F)
#' }

fsd.spca.filters = function (F, Npc = 1, L = 3)
{
  if (!inherits(F, "fsd.freqdom"))
    stop("F needs to be a spectral density operator")

  L = floor(L)

  if (length(L) == 1)
    L = rep(L, length(F$freq[[1]]))

  if (length(L) != length(F$freq[[1]]))
    stop("Dimensions incompatible")

  dims = dim(F$operators[[1]])[1]
  Npc = min(Npc, dims)

  # Change of base
  if (!is.null(F$basis)) {
    B = inprod(F$basis, F$basis)
  } else {
    B = diag(dim(F$operators[[1]])[1])
  }
  B.root = eigen(B)$vectors %*% diag(sqrt(eigen(B)$values)) %*%
           t(eigen(B)$vectors)
  B.root.minus = solve(B.root)
  for (i in 1:length(F$operators)) {
    F$operators[[i]] = B.root %*% F$operators[[i]] %*% B.root
  }

  # Compute eigenvector system
  Esys = lapply(F$operators, eigen, symmetric = TRUE)
  Evecs = sapply(Esys, function(w) {w$vectors[ ,1:Npc, drop=FALSE]},
             simplify="array")

  # Compute portion of variance explained
  Evals = sapply(Esys, '[[', "values", simplify="array")
  # Integrate
  Evals.int = apply(pmax(Re(Evals), 0), 1, sum)
  names(Evals.int) = paste("PC", 1:dim(Evals)[1])

  sfpca.var = Evals.int/sum(Evals.int)

  # Align the eigenvectors in the direction of v
  v = rep(1, dims)
  phiv = apply(Evecs, 2:3, "%*%", v)
  phiv[abs(phiv) < 1e-14] = 1
  phiv = phiv/abs(phiv)

  Evecs = Evecs / rep(phiv, each = dims)

  # Compute the Fourier inverse
  laglist = unfold(lapply(L, function(h) {-h:h}))

  filters = fsd.fourier.inverse(Evecs, F$freq, laglist)

  # Change of base
  for (i in 1:length(filters$operators)) {
    filters$operators[[i]] = B.root.minus %*% filters$operators[[i]]
  }

  filters$basis = F$basis

  filters$var = sfpca.var

  # Check how much the PC filters cover
  phinorms = sapply(filters$operators, function(M) {diag(t(M) %*% B %*% M)})

  if (Npc == 1)
    phinorms = sum(phinorms)
  else
    phinorms = apply(phinorms, 1, sum)

  if (any(phinorms < 0.9))
    warning("The chosen clipping parameter L is low.\n",
            "Summed squared norms for the PC filters are only ",
            toString(round(phinorms, 3)))

  filters
}
