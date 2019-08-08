#' Creating a functional spatial filter
#'
#' This function creates a functional spatial filter
#' @param laglist a list of lags.
#' @param operatorlist a list of operators, i.e. matrices.
#' @param basisobj a functional basis object defining the basis.
#' @return the \code{fsd.filter} object.
#' @keywords fsd
#' @export
#' @examples
#' \dontrun{
#' fsd.filter(laglist, operatorlist)
#' }

fsd.filter = function (laglist, operatorlist, basisobj = NULL)
{
  if (length(laglist) != length(operatorlist))
    stop("There need to be as many operator matrices as lags.")

  A = list(laglist = laglist,
           operators = lapply(operatorlist, as.matrix),
           basis = basisobj)
  class(A) = "fsd.filter"
  A
}
