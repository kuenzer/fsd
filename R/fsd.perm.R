#' Permute the dimensions of a Functional Spatial Data Grid
#'
#' This function is used to permute the order of the dimensions of the spatial
#' grid of some functional spatial data.
#' @param fsdobj the functional spatial data.
#' @param perm the subscript permutation vector, i.e. a permutation of the
#'   integers 1:r.
#' @seealso \code{\link{fsd.fd}}
#' @keywords fsd
#' @export
#' @examples
#' data("temp")
#' fsd.plot.data(fsd.perm(temp, c(2, 1, 3)))

fsd.perm = function (fsdobj, perm = NULL)
{
  if (!inherits(fsdobj, "fsd.fd"))
    stop("fsdobj needs to be an fsd object")

  if (length(dim(fsdobj$coefs)) == 3)
    perm = 2:1
  else if (length(dim(fsdobj$coefs)) <= 2)
    return(fsdobj)

  if (!is.numeric(perm) || length(perm) != length(dim(fsdobj$coefs)) - 1)
    stop("Invalid permutation")

  perm = c(1, 1 + perm)

  fsdobj$coefs = aperm(fsdobj$coefs, perm)
  fsdobj$fdnames = fsdobj$fdnames[perm]

  fsdobj
}
