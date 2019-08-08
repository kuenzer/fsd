#' Check if functional spatial filter is unilateral
#'
#' This function is used to verify if a functional spatial filter is unilateral.
#' @param A the functional spatial filter.
#' @param exclude.origin whether to exclude the origin from the definition of
#'   unilaterality.
#' @return a Boolean indicating whether A is unilateral or not.
#' @seealso \code{\link{fsd.filter}}
#' @keywords fsd
#' @export
#' @examples
#' \dontrun{
#' fsd.filter.is.unilateral(A)
#' }

fsd.filter.is.unilateral = function (A, exclude.origin = FALSE)
{
  if (!inherits(A, "fsd.filter"))
    stop("Input needs to be a filter object.")

  r = length(A$laglist[[1]])
  lags = matrix(unlist(A$laglist), nrow = r)

  if (exclude.origin)
    if (any(apply(lags == 0, 2, all)))
      return(FALSE)

  for (k in 1:ncol(lags)) {
    for (i in 1:r) {
      if (lags[i,k] < 0)
        return(FALSE)
      else if (lags[i,k] > 0)
        break
    }
  }

  return(TRUE)
}
