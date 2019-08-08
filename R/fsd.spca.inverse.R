#' Reconstruct the original functional data from the filters and the score
#'
#' This function is used to reconstruct the original functional data from the
#' filters and the score. It applies the transposed matrices of the filter
#' \code{A} to the scores.
#' @param A the filters.
#' @param scores the scores as an array.
#' @param mean.X the mean of X.
#' @return X the reconstructed functional data
#' @seealso \code{\link{fsd.spca}}, \code{\link{fsd.spca.scores}}
#' @keywords fsd
#' @export
#' @examples
#' \dontrun{
#' fsd.spca.inverse(A, scores)
#' }

fsd.spca.inverse = function (A, scores, mean.X = NULL)
{
  if (!inherits(A, "fsd.filter"))
    warning("A is not of the class fsd.filter")

  if (dim(scores)[1] > dim(A$operators[[1]])[2])
    stop("Dimensions are not compatible: Number of PC too small")
  if (length(dim(scores)) - 1 != length(A$laglist[[1]]))
    stop("Dimensions are not compatible: Filter lags and scores do not fit")

  if (dim(scores)[1] < dim(A$operators[[1]])[2]) {
    warning("Less scores than PC filters provided. ",
            "Remaining scores are assumed to be zero.")
    dimdiff = dim(A$operators[[1]])[2] - dim(scores)[1]
    scores = apply(scores, 2:length(dim(scores)),
                   function(y) {c(y, rep(0, dimdiff))})
  }

  A2 = A
  A2$laglist = lapply(A2$laglist, "-")
  A2$operators = lapply(A2$operators, t)
  A2$basis = NULL

  X = fsd.spca.scores(scores, A2)

  if (!is.null(mean.X)) {
    if (inherits(mean.X, "fd"))
      mean.X = mean.X$coefs
    if (length(mean.X) == dim(X)[1]) {
      X = X + as.vector(mean.X)
    } else {
      warning("Mean X does not have the correct dimension")
    }
  }

  dimnames(X)[[1]] = sprintf("f%s", 1:dim(X)[1])

  if (!is.null(A$basis)) {
    dimnames(X)[[1]] = A$basis$names
    X = fsd.fd(X, A$basis)
  }

  X
}
