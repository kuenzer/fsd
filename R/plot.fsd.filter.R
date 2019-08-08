#' Plot functional spatial filters
#'
#' This function is used to plot functional spatial filters,  along a
#' two-dimensional plane of indices.
#' @param x the filter.
#' @param gridsize the number of grid points on which to evaluate the functional
#'   data.
#' @param Npc number of filters to plot.
#' @param Lmax maximum lag to plot.
#' @param set.dims the coordinates to be fixed. A vector of dimension r. NaN
#'   represents not fixed components.
#' @param ... graphical parameters.
#' @seealso \code{\link{fsd.plot.filters}}
#' @keywords fsd
#' @export
#' @examples
#' \dontrun{
#' plot.fsd.filter(x)
#' }

plot.fsd.filter = function (x, ..., gridsize = NULL, Npc = 3, Lmax = 3,
                            set.dims = NULL)
{
  arglist <- list(...)

  fsd.plot.filters(x, gridsize = gridsize, Npc = Npc, Lmax = Lmax,
                   set.dims = set.dims, main = arglist$main,
                   xlab = arglist$xlab, ylab = arglist$ylab)
}
