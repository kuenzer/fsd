#' Compute the norm of functional data
#'
#' This function is used to compute the norm of functional data.
#' @param X the functional spatial data.
#' @return an array of the norms.
#' @keywords fsd
#' @export
#' @examples
#' \dontrun{
#' fsd.norm(X)
#' }

fsd.norm = function (X)
{
  if (!inherits(X, "fd"))
    warning("X is not a functional data object")

  B = inprod(X$basis, X$basis)
  r = length(dim(X$coefs)) - 1

  sqrt(apply(X$coefs, 1+1:r, function(x) {t(x) %*% B %*% x}))
}
