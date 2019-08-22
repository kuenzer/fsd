#' Compute the Fourier transform
#'
#' This function is used to calculate the Fourier transform
#' @param covs the covariances to be transformed
#' @param laglist list of lags.
#' @param thlist list of frequencies.
#' @return the Fourier transform of the covariances.
#' @keywords fsd
#' @importFrom Rcpp evalCpp
#' @examples
#' \dontrun{
#' fsd.fourier(covs, thlist, lags)
#' }

fsd.fourier = function (covs, laglist, thlist)
{
  if (length(thlist[[1]]) != length(laglist[[1]]))
    stop("Lags and frequencies need to have the same dimension.")
  if (dim(covs)[3] != length(laglist))
    stop("There need to be as many lags as matrices of covariances")

  A = fsd_rcpp_fourier(covs,
                       matrix(unlist(laglist),
                              nrow = length(laglist[[1]])),
                       matrix(unlist(thlist),
                              nrow = length(thlist[[1]])))
  dim(A) = c(dim(covs)[1:2], length(thlist))

  # A needs to be transformed from an array into a list of matrices
  A = sapply(1:dim(A)[3], function(i){A[,,i]}, simplify = FALSE)

  Aobj = list(operators = A, thlist = thlist)
  Aobj
}
