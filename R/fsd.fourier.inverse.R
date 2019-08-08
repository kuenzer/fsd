#' Compute the Fourier inverse
#'
#' This function is used to calculate the Fourier inverse.
#' @param evecs the eigenvectors to be expanded.
#' @param thlist list of frequencies.
#' @param laglist list of lags.
#' @return the filter obtained by applying the Fourier inverse on evecs.
#' @keywords fsd
#' @examples
#' \dontrun{
#' fsd.fourier.inverse(evecs, thlist, lags)
#' }

fsd.fourier.inverse = function (evecs, thlist, laglist)
{
  if (length(thlist[[1]]) != length(laglist[[1]]))
    stop("Lags and frequencies need to have the same dimension.")
  if (dim(evecs)[3] != length(thlist))
    stop("There need to be as many frequencies as matrices of eigenvectors.")

  A = Re(fsd_rcpp_fourier_inverse(evecs,
                                  matrix(unlist(thlist),
                                         nrow = length(thlist[[1]])),
                                  matrix(unlist(laglist),
                                         nrow = length(laglist[[1]]))))
  dim(A) = c(dim(evecs)[1:2], length(laglist))

  # A needs to be transformed from an array into a list of matrices
  A = sapply(1:dim(A)[3], function(i){A[,,i]}, simplify = FALSE)

  Aobj = list(operators = A, laglist = laglist)
  class(Aobj) = "fsd.filter"
  Aobj
}
