#' Plot functional spatial data
#'
#' This function is used to plot functional spatial data on a grid of dimension
#' r,  along a two-dimensional plane of indices.
#' @param X the functional spatial data. Either an fd object or an array, or a
#'   list thereof.
#' @param basisobj the basis of the functional data.
#' @param gridsize the number of grid points on which to evaluate the functional
#'   data.
#' @param set.dims the coordinates to be fixed. A vector of dimension r. NaN
#'   represents not fixed components.
#' @param xlim limits for the coordinates on the x-axis.
#' @param ylim limits for the coordinates on the y-axis.
#' @param main title of the plot.
#' @param xlab label for the x-axis.
#' @param ylab label for the y-axis.
#' @seealso \code{\link{plot.fsd.fd}}, \code{\link{fsd.plot.filters}},
#'   \code{\link{fsd.plot.covariance}}
#' @keywords fsd
#' @import stats
#' @export
#' @examples
#' data("temp")
#' # plot the data for the third year in the dataset (i.e. 1981)
#' fsd.plot.data(temp, set.dims = c(NaN, NaN, 3))

fsd.plot.data = function (X, basisobj = NULL, gridsize = NULL,
                          set.dims = NULL, xlim = NULL, ylim = NULL,
                          main = NULL, xlab = NULL, ylab = NULL)
{
  if (!is.list(X) || inherits(X, "fd"))
    X = list(X)

  fdnames = NULL
  basislist = lapply(X, function(Z) {if (inherits(Z, "fd")) Z$basis else NULL})
  namelist = lapply(X, function(Z) {if (inherits(Z, "fd")) Z$fdnames else NULL})
  X = lapply(X, function(Z) {if (inherits(Z, "fd")) Z$coefs else Z})

  if (length(X) > 1) {
    if (sd(sapply(lapply(X, dim), length)) > 0)
      stop("Inconsistent dimensions: all X need to have the same dimensions")
    else if (any(apply(sapply(X, dim), 1, sd) > 0))
      stop("Inconsistent dimensions: all X need to have the same dimensions")
  }

  isfsdobj = !sapply(basislist,is.null)
  if (any(isfsdobj)) {
    basisobj = basislist[[which(isfsdobj)[1]]]
    fdnames = namelist[[which(isfsdobj)[1]]]
  }
  for (k in which(isfsdobj)) {
    if (!(basisobj == basislist[[k]]))
      stop("Inconsistent bases: all X need to have the same basis (or none)")
    if (any(names(fdnames) != names(namelist[[k]])))
      warning("Inconsistent dimension naming: ",
              "not all X have the same dimension names")
  }

  if (is.null(basisobj))
    basisobj = create.bspline.basis(nbasis = dim(X[[1]])[1], norder = 1)

  r = length(dim(X[[1]])) - 1

  if (is.null(set.dims)) {
    if (r == 1)
      set.dims = c(NaN)
    else
      set.dims = c(NaN, NaN)
    if (r > 2) {
      set.dims = c(set.dims, rep(1, r - 2))
    }
  }
  if (length(set.dims) != r)
    stop("The parameter set.dims needs the same length as the dimensions of X")

  if (sum(is.nan(set.dims)) > 2)
    stop("You can only plot two dimensions")

  freedims = which(is.nan(set.dims))
  if (length(freedims) == 0)
    stop("Please provide a dimension to plot")
  only.1d = (length(freedims) == 1)

  if (is.null(xlim))
    xlim = c(1, dim(X[[1]])[1+freedims[1]])
  if (is.null(ylim)) {
    if (!only.1d)
      ylim = c(1, dim(X[[1]])[1+freedims[2]])
    else
      ylim = c(1, 1)
  }

  xlim = floor(xlim)
  ylim = floor(ylim)

  xlim[1] = max(xlim[1], 1)
  xlim[2] = min(xlim[2], dim(X[[1]])[1+freedims[1]])

  if (!only.1d) {
    ylim[1] = max(ylim[1], 1)
    ylim[2] = min(ylim[2], dim(X[[1]])[1+freedims[2]])
  } else {
    ylim = c(1, 1)
  }

  inds = as.list(set.dims)
  inds[[freedims[1]]] = xlim[1]:xlim[2]
  if (!only.1d)
    inds[[freedims[2]]] = ylim[1]:ylim[2]

  Xplot = lapply(X, function(Z) {Y = do.call(what = apply,
                                           args = c(list(X = Z, MARGIN = 1,
                                                         FUN = '['),
                                                    inds))
                                 if (!is.vector(Y))
                                   Y = t(Y)
                                 return(Y)})

  if (only.1d)
    positions = as.list(xlim[1]:xlim[2])
  else
    positions = unfold(list(xlim[1]:xlim[2], ylim[1]:ylim[2]))

  if (is.null(main)) {
    main = "Data Plot"

    if (any(!is.nan(set.dims))) {
      main = c(main,
               paste("with indices (",
                     paste(gsub(NaN, "*", set.dims), collapse = " , "), ")"))
    }
  }

  if (is.null(xlab)) {
    if (!is.null(fdnames)) {
      xlab = names(fdnames)[1+freedims[1]]
    } else {
      xlab = paste("Dimension", freedims[1])
    }
  }

  if (is.null(ylab)) {
    ylab = ""
    if (!only.1d) {
      if (!is.null(fdnames)) {
        ylab = names(fdnames)[1+freedims[2]]
      } else {
        ylab = paste("Dimension", freedims[2])
      }
    }
  }

  if (is.null(gridsize))
    gridsize = max(50, ceiling(c(500, 10 * basisobj$nbasis)/(diff(xlim) + 1)))

  fsd.z.plot(Zlist = Xplot, positions = positions,
             basisobj = basisobj, gridsize = gridsize, main = main,
             xlab = xlab, ylab = ylab)
}
