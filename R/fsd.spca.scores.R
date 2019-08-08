#' Compute the Scores for Spatial Functional Data
#'
#' This function computes the scores. It applies the matrices of the filter
#' \code{A} to the data \code{X}.
#' @param X the functional spatial data. Either an fd object or an array.
#' @param A the filter used to compute the scores.
#' @param na.ignore whether missing values in the data X should be ignored.
#' @keywords fsd
#' @return the array of the scores.
#' @seealso \code{\link{fsd.spca.inverse}}
#' @export
#' @examples
#' \dontrun{
#' fsd.spca.scores(X, A)
#' }

fsd.spca.scores = function (X, A, na.ignore = FALSE)
{
  if (inherits(X, "fd"))
    X = X$coefs

  r = length(dim(X)) - 1
  d = dim(X)[1]

  if (dim(X)[1] != dim(A$operators[[1]])[1])
    stop("Dimensions are not compatible: Number of basis functions varies")
  if (r != length(A$laglist[[1]]))
    stop("Dimensions are not compatible: Filter lags and data do not fit")

  # Change of base
  if (!is.null(A$basis)) {
    B = inprod(A$basis, A$basis)
    X = apply(X, 1+1:r, `%*%`, B)
  }

  Y = array(0, dim = c(dim(A$operators[[1]])[2], dim(X)[1+1:r]),
            dimnames = c(list(paste("PC", 1:dim(A$operators[[1]])[2])),
                         dimnames(X)[1+1:r]))

  if (na.ignore) {
    Xna = apply(is.na(X), 1+1:r, any)
    Y[1:dim(Y)[1] + dim(Y)[1] * rep(which(Xna) - 1, each = dim(Y)[1])] = NA
    X[1:dim(X)[1] + dim(X)[1] * rep(which(Xna) - 1, each = dim(X)[1])] = 0
  }

  for (i in 1:length(A$operators)) {
    if (any(dim(X)[1+1:r] <= abs(A$laglist[[i]])))
      next

    inds1 = as.list(rep(TRUE, 1+r))
    inds2 = inds1

    for (j in 1:length(A$laglist[[i]])) {
      if (A$laglist[[i]][j] > 0) {
        inds1[[1+j]] = -(dim(X)[1+j] + 1 - (1:abs(A$laglist[[i]][j])))
        inds2[[1+j]] = -(1:abs(A$laglist[[i]][j]))
      }
      else if (A$laglist[[i]][j] < 0) {
        inds1[[1+j]] = -(1:abs(A$laglist[[i]][j]))
        inds2[[1+j]] = -(dim(X)[1+j] + 1 - (1:abs(A$laglist[[i]][j])))
      }
    }

    Y = do.call(what = `[<-`,
                args = c(list(Y), inds2,
                         list(do.call(what = `[`,
                                      args = c(list(Y), inds2, drop = FALSE)) +
                              as.vector(t(A$operators[[i]]) %*%
                                        matrix(do.call(what = `[`,
                                                       args = c(list(X), inds1,
                                                                drop = FALSE)),
                                               nrow = d)))))
  }

  Y
}
