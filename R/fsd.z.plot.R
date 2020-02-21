#' Plot functional data on a grid
#'
#' This function plots functional data on a grid.
#' @param Zlist list of arrays with coefficients.
#' @param positions the positions of the single curves in space.
#' @param basisobj the basis of the functional data.
#' @param gridsize the number of grid points on which to evaluate the curves.
#' @param ... other arguments to pass to the plot function
#' @keywords fsd
#' @import stats fda
#' @importFrom graphics plot lines abline axis
#' @examples
#' \dontrun{
#' fsd.z.plot(Zlist, positions, basisobj)
#' }

fsd.z.plot = function (Zlist, positions, basisobj = NULL, gridsize = 50, ...)
{
  if (is.null(basisobj))
    basisobj = create.bspline.basis(nbasis = dim(Zlist[[1]])[1], norder = 1)

  egr = 0:gridsize / gridsize
  egr2 = egr * diff(basisobj$rangeval) + basisobj$rangeval[1]

  Zev = sapply(Zlist, function(x) {eval.fd(egr2, fd(coef = x,
                                                    basisobj = basisobj))},
               simplify = "array")

  if (is.vector(Zev))
    Zev = matrix(Zev, ncol = 1)

  Zev = Zev/max(abs(Zev), na.rm = TRUE) /2 *0.95

  only.1d = (length(positions[[1]]) == 1)

  lags1 = unique(sapply(positions, '[', 1))
  lags2 = 0
  if (!only.1d)
    lags2 = unique(sapply(positions, '[', 2))

  xlim = c( min(lags1) , max(lags1) )
  ylim = c( min(lags2) , max(lags2) )

  plot(NULL, NULL,
       xlim = xlim + c(-0.5, 0.5), ylim = ylim + c(-0.5, 0.5),
       xaxt = 'n', yaxt = 'n', frame.plot = FALSE, ...)

  d1p = pretty(xlim[1]:xlim[2], n = 10)
  axis(1, at = d1p, labels = d1p, tick = FALSE)
  if (!only.1d) {
    d2p = pretty(ylim[1]:ylim[2], n = 10)
    axis(2, at = d2p, labels = d2p, tick = FALSE)
  }

  abline(h = lags2, v = lags1, col = 8, lty = 3)

  for (m in length(Zlist):1) {
    for (i in 1:dim(Zev)[2]) {
      lines(positions[[i]][1] + (egr -0.5)*0.8,
            ifelse(!only.1d, positions[[i]][2], 0) + Zev[ , i, m], col = m)
    }
  }
}
