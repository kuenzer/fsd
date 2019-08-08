#' Plot functional spatial data
#'
#' This function is used to plot functional spatial data on a grid of dimension
#' r,  along a two-dimensional plane of indices.
#' @param x the functional spatial data.
#' @param gridsize the number of grid points on which to evaluate the functional
#'   data.
#' @param set.dims the coordinates to be fixed. A vector of dimension r. NaN
#'   represents not fixed components.
#' @param ... graphical parameters.
#' @seealso \code{\link{fsd.plot.data}}
#' @keywords fsd
#' @export
#' @examples
#' data("temp")
#' plot(temp)

plot.fsd.fd = function (x, ..., gridsize = NULL, set.dims = NULL)
{
  arglist <- list(...)

  fsd.plot.data(x, gridsize = gridsize, set.dims = set.dims,
                xlim = arglist$xlim, ylim = arglist$ylim, main = arglist$main,
                xlab = arglist$xlab, ylab = arglist$ylab)
}
