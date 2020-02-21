#' Plot the autocovariance operator of functional spatial data
#'
#' This function is used to plot the autocovariance operator of functional
#' spatial data along a two-dimensional plane of lags.
#' @param X the functional spatial data. Either an fd object or an array.
#' @param basisobj the basis of the functional data.
#' @param q the maximal lag to be plotted.
#' @param gridsize the resolution of the grid on which the kernel will be
#'   evaluated.
#' @param set.dims the coordinates to be fixed. A vector of dimension r. NaN
#'   represents not fixed components.
#' @param plot.cor a boolean indicating whether the correlation or the
#'   covariance should be plotted.
#' @param na.rm whether or not missing values should be ignored.
#' @param main title of the plot.
#' @param xlab label for the x-axis.
#' @param ylab label for the y-axis.
#' @param ... other arguments to pass to the plot function
#' @seealso \code{\link{fsd.plot.data}}
#' @keywords fsd
#' @import fda
#' @importFrom graphics filled.contour abline axis
#' @export
#' @examples
#' data("temp")
#' fsd.plot.covariance(temp, q = 3, plot.cor = TRUE, zlim = c(-1,1))

fsd.plot.covariance = function (X, basisobj = NULL, q = 0, gridsize = NULL,
                                set.dims = NULL, plot.cor = FALSE,
                                na.rm = FALSE, main = NULL,
                                xlab = NULL, ylab = NULL, ...)
{
  fdnames = NULL
  if (inherits(X, "fd")) {
    basisobj = X$basis
    fdnames = X$fdnames
    X = X$coefs
  }

  if (is.null(basisobj)) {
    basisobj = create.bspline.basis(nbasis = dim(X)[1], norder = 1)
  }

  r = length(dim(X)) - 1

  if (is.null(set.dims)) {
    if (r == 1)
      set.dims = c(NaN)
    else
      set.dims = c(NaN, NaN)
    if (r > 2) {
      set.dims = c(set.dims, rep(0, r - 2))
    }
  }
  if (length(set.dims) != r)
    stop("Need right dimensions")

  if (sum(is.nan(set.dims)) > 2)
    stop("You can only plot two dimensions")

  freedims = which(is.nan(set.dims))
  if (length(freedims) == 0)
    stop("Please provide a dimension to plot")
  only.1d = (length(freedims) == 1)

  q = floor(q)
  if (length(q) == 1) {
    q = rep(q, 2)
  } else if (length(q) > 2) {
    q = q[1:2]
  }

  q[1] = min(q[1], dim(X)[1+freedims[1]]-1)
  if (!only.1d)
    q[2] = min(q[2], dim(X)[1+freedims[2]]-1)
  else
    q[2] = 0

  q0 = rep(0, r)
  q0[freedims] = q[1:length(freedims)]

  inds = lapply(q0, function(w) {-w:w})
  inds[-freedims] = set.dims[!is.nan(set.dims)]
  hlist = unfold(inds)

  singlebasis = (basisobj$nbasis == 1)

  if (singlebasis)
    gridsize = 2

  if (is.null(gridsize))
    gridsize = max(40, ceiling(c(200, 8 * basisobj$nbasis)/(2*min(q) + 1)))

  X = X - apply(X, 1, mean, na.rm = na.rm)

  Clist = lapply(hlist, FUN = fsd.covariance, Y = X,
                 centered = TRUE, na.rm = na.rm)

  gr = (1:gridsize - 0.5) / gridsize
  egr = gr * diff(basisobj$rangeval) + basisobj$rangeval[1]
  lgr = length(gr)

  Cev = lapply(Clist,
               function(C) {eval.bifd(egr, egr,
                                      bifd(coef = C, sbasisobj = basisobj,
                                           tbasisobj = basisobj))})

  if (plot.cor) {
    zeroindex = 1
    if (length(hlist) > 1)
      zeroindex = which(sapply(hlist, function(h) {all(h == 0)}))
    if (length(zeroindex) == 0) {
      C0 = fsd.covariance(Y = X, h = rep(0, r), centered = TRUE)
      vars = diag( eval.bifd(egr, egr, bifd(coef = C0, sbasisobj = basisobj,
                                            tbasisobj = basisobj)) )
    } else {
      vars = diag(Cev[[zeroindex]])
    }
    nm = 1/sqrt(vars%*%t(vars))
    Cev = lapply(Cev, "*", nm)
  }
  bigC = array(dim = lgr* (2*q+1))

  for (k in 1:length(Clist)) {
    bigC[(q[1] + hlist[[k]][freedims[1]]) * lgr + 1:lgr,
         ifelse(!only.1d, q[2] + hlist[[k]][freedims[2]],
                0) * lgr + 1:lgr] = Cev[[k]]
  }
  if (plot.cor)
    bigC = pmin(pmax(bigC, -1), 1)


  if (is.null(main)) {
    if (!plot.cor)
      main = "Covariance Operator Plot"
    else
      main = "Correlation Operator Plot"

    if (any(!is.nan(set.dims))) {
      main = c(main,
               paste("with lags (",
                     paste(gsub(NaN, "*", set.dims), collapse = " , "), ")"))
    }
  }

  if (is.null(xlab)) {
    if (!is.null(fdnames)) {
      xlab = paste("Lag for", names(fdnames)[1+freedims[1]])
    } else {
      xlab = paste("Lag in dimension", freedims[1])
    }
  }

  if (is.null(ylab)) {
    ylab = ""
    if (!only.1d) {
      if (!is.null(fdnames)) {
        ylab = paste("Lag for", names(fdnames)[1+freedims[2]])
      } else {
        ylab = paste("Lag in dimension", freedims[2])
      }
    }
  }

  filled.contour(rep(-q[1]:q[1], each=lgr) + gr, rep(-q[2]:q[2], each=lgr) + gr,
                 bigC, xlim = c(-q[1],q[1]+1), ylim = c(-q[2],q[2]+1),
                 asp = 1, xlab = xlab, ylab = ylab,
                 plot.axes = {d1p = -q[1]:q[1]
                 if (length(d1p) > 10)
                   d1p = pretty(d1p, n = 10)
                 axis(1, at = d1p+0.5, labels = d1p, tick = singlebasis)
                 if (!only.1d) {
                   d2p = -q[2]:q[2]
                   if (length(d2p) > 10)
                     d2p = pretty(d2p, n = 10)
                   axis(2, at = d2p+0.5, labels = d2p, tick = singlebasis)
                 }
                 if (!singlebasis)
                   abline(v = -q[1]:(q[1]+1),
                          h = -q[2]:(q[2]+1), col = 8)
                 },
                 frame.plot = FALSE, main = main, ...)
}
