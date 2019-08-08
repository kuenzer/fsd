#' Plot functional spatial filters
#'
#' This function is used to plot functional spatial filters,  along a
#' two-dimensional plane of indices.
#' @param A the filters.
#' @param basisobj the basis of the functional data.
#' @param gridsize the number of grid points on which to evaluate the functional
#'   data.
#' @param Npc number of filters to plot.
#' @param Lmax maximum lag to plot.
#' @param set.dims the coordinates to be fixed. A vector of dimension r. NaN
#'   represents not fixed components.
#' @param main title of the plot.
#' @param xlab label for the x-axis.
#' @param ylab label for the y-axis.
#' @seealso \code{\link{plot.fsd.filter}}, \code{\link{fsd.plot.data}},
#'   \code{\link{fsd.plot.covariance}}
#' @keywords fsd
#' @export
#' @examples
#' \dontrun{
#' fsd.plot.filters(A)
#' }

fsd.plot.filters = function (A, basisobj = NULL, gridsize = NULL,
                             Npc = 3, Lmax = 3, set.dims = NULL,
                             main = NULL, xlab = NULL, ylab = NULL)
{
  if (inherits(A, "fsd.filter")) {
    basisobj = A$basis
  }

  if (is.null(basisobj))
    basisobj = create.bspline.basis(nbasis = dim(A$operators[[1]])[1],
                                    norder = 1)

  if (is.null(set.dims)) {
    if (length(A$laglist[[1]]) == 1)
      set.dims = c(NaN)
    else
      set.dims = c(NaN, NaN)
    if (length(A$laglist[[1]]) > 2) {
      set.dims = c(set.dims, rep(0, length(A$laglist[[1]]) - 2))
    }
  }
  if (length(set.dims) != length(A$laglist[[1]]))
    stop("Need right dimensions")

  if (sum(is.nan(set.dims)) > 2)
    stop("You can only plot two dimensions")

  freedims = which(is.nan(set.dims))
  if (length(freedims) == 0)
    stop("Please provide a dimension to plot")
  only.1d = (length(freedims) == 1)

  Npc = min(Npc, dim(A$operators[[1]])[2])

  lagrange = apply(matrix(unlist(A$laglist), nrow = length(A$laglist),
                          byrow = TRUE), 2, range)

  if (any(set.dims > lagrange[2,], na.rm = TRUE) ||
      any(set.dims < lagrange[1,], na.rm = TRUE))
    stop("Please choose an available lag")

  plotl = sapply(A$laglist,
                 function(w) {(max(abs((is.nan(set.dims)) * w)) <= Lmax) &&
                     all(w[!is.nan(set.dims)] == set.dims[!is.nan(set.dims)])})
  laglist = lapply(A$laglist[plotl], '[', freedims)

  Alist = lapply(seq(Npc), function(m) {sapply(A$operators[plotl],
                                               function(M) {M[,m]})})

  if (is.null(main)) {
    main = "SFPC Filter Plot"

    if (any(!is.nan(set.dims))) {
      main = c(main,
               paste("with lags (",
                     paste(gsub(NaN, "*", set.dims), collapse = " , "), ")"))
    }
  }

  if (is.null(xlab)) {
    xlab = paste("Lags in dimension", freedims[1])
  }

  if (is.null(ylab)) {
    ylab = ""
    if (!only.1d) {
      ylab = paste("Lags in dimension", freedims[2])
    }
  }

  if (is.null(gridsize))
    gridsize = max(50, ceiling(c(500, 10 * basisobj$nbasis) /
                                 (diff(range(sapply(laglist, '[', 1))) + 1)))

  fsd.z.plot(Zlist = Alist, positions = laglist,
             basisobj = basisobj, gridsize = gridsize, main = main,
             xlab = xlab, ylab = ylab)
}
