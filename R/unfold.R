#' Unfold indices into a high-dimensional grid
#'
#' This function is used to expand indices into a grid. Works similar to
#' \code{\link[base]{expand.grid}}.
#' @param indexlist indices to be expanded. Either a list with each member
#'   representing a dimension, or a vector.
#' @param r the dimension of the grid to be created. This is only needed if
#'   indexlist is a vector.
#' @param return.array a boolean indicating whether an array should be returned
#'   or a list.
#' @return a list or an array with the coordinates of the grid.
#' @keywords fsd
#' @examples
#' \dontrun{
#' fsd:::unfold(list(11:12, 21:23), return.array = TRUE)
#' }

unfold = function (indexlist, r = NULL, return.array = FALSE)
{
  if (!is.list(indexlist)) {
    tmpl = list()
    for (i in 1:r) {
      tmpl[[i]] = indexlist
    }
    indexlist = tmpl
  }

  r = length(indexlist)

  gridlist = as.list(indexlist[[1]])

  if (r > 1) {
    for (i in 2:r) {
      gridlist = unlist(lapply(indexlist[[i]],
                               function(x) {lapply(gridlist,
                                                   function(y) {c(y, x)})}),
                        recursive=FALSE)
    }
  }

  if (return.array) {
    return( matrix(unlist(gridlist), nrow = length(gridlist), byrow = TRUE) )
  } else {
    return(gridlist)
  }
}
