#' Simulate a spatial functional ARMA process
#'
#' This function is used to simulate a sample from a spatial functional ARMA
#' (SFARMA) process using a normal or t-distributed noise process.
#'
#' Please note that convergence in the fixed point iteration depends strongly on
#' the sample size and the structure of the AR filter. For certain choices of
#' parameters, even more than \code{n} iterations may be needed.
#' @param n the sample size.
#' @param Sigma the covariance matrix of the noise.
#' @param ARfilter the functional spatial filter of the autoregressive part.
#' @param MAfilter the functional spatial filter of the moving-average part.
#' @param burnin the number of grid points to add on each side of the sample to
#'   achieve stationarity.
#' @param basisobj the basis of the functional data.
#' @param noise the noise distribution to use. Either "normal" or an integer for
#'   the degrees of freedom of a multivariate Student's t-distribution.
#' @param do.fixed.point whether to always use a fixed point iteration.
#' @param max.iter the maximum number of iterations in the fixed point
#'   iteration.
#' @param eps the stopping criterion for the fixed point iteration.
#' @return one sample of the SFARMA process.
#' @seealso \code{\link{fsd.fd}}
#' @keywords fsd
#' @import fda mvtnorm
#' @export
#' @examples
#' \dontrun{
#' fsd.sfarma(n, Sigma, ARfilter, MAfilter, basisobj)
#' }

fsd.sfarma = function (n, Sigma = NULL, ARfilter = NULL, MAfilter = NULL,
                       burnin = 30, basisobj = NULL, noise = "normal",
                       do.fixed.point = FALSE, max.iter = 50, eps = 1e-5)
{

  if (is.null(Sigma)) {
    if (!is.null(ARfilter))
      d = dim(ARfilter$operators[[1]])[1]
    else if (!is.null(MAfilter))
      d = dim(MAfilter$operators[[1]])[1]
    else
      stop("Not enough parameters")
    Sigma = diag(exp(- (1:d -1) / 10))
  }

  Sigma = as.matrix(Sigma)

  d = dim(Sigma)[1]
  r = length(n)

  if (is.null(basisobj)) {
    basisobj = create.bspline.basis(nbasis = d, norder = 1)
  }

  if (!inherits(basisobj, "basisfd"))
    stop("basisobj needs to be a basis.")

  if (!is.null(ARfilter) && !inherits(ARfilter, "fsd.filter"))
    stop("ARfilter needs to be a filter object.")
  if (!is.null(MAfilter) && !inherits(MAfilter, "fsd.filter"))
    stop("MAfilter needs to be a filter object.")

  if (!is.numeric(noise) || noise == 0)
    Z = array(t(rmvnorm(prod(n + 2*burnin), sigma = Sigma)),
              dim = c(d, n + 2*burnin))
  else
    Z = array(t(rmvt(prod(n + 2*burnin), sigma = Sigma, df = noise)),
              dim = c(d, n + 2*burnin))

  if (!is.null(MAfilter))
    Z = fsd.spca.scores(X = Z, A = MAfilter)

  X = Z

  if (!is.null(ARfilter)) {
    if (!fsd.filter.is.unilateral(ARfilter) && !do.fixed.point) {
      message("Multilateral AR-filter necessitates fixed point iterations...")
      do.fixed.point = TRUE
    }
    if (!do.fixed.point) {
      message("Unilateral AR-filter permits ordinary recursive computation.")
      maxARlag = sapply(split(unlist(ARfilter$laglist),f = seq(r)), max)
      for (i in 1:prod(n + 2*burnin)) {
        ind = arrayInd(i, n + 2*burnin)
        if (any(ind < 1 + maxARlag))
          next

        for (k in 1:length(ARfilter$operators)) {
          if (all(ARfilter$laglist[[k]] == 0))
            next
          X = do.call(what = `[<-`,
                      args = c(list(X), TRUE, ind,
                               list(do.call(what = `[`, args = c(list(X), TRUE,
                                                                 ind)) +
                                      ARfilter$operators[[k]] %*%
                                      do.call(what = `[`,
                                              args = c(list(X), TRUE,
                                                       ind - ARfilter$laglist[[k]])))))
        }
      }
    } else {
      message("Starting fixed point iterations with eps = ", eps,
              " and max.iter = ", max.iter)
      convergent = FALSE
      for (i in 1:max.iter) {
        Xold = X
        X = fsd.spca.scores(X = X, A = ARfilter) + Z
        difference = max(fsd.norm(fsd.fd(X - Xold, basisobj = basisobj)))
        if (difference < eps) {
          convergent = TRUE
          break
        }
      }
      if (convergent)
        message("Convergence reached after ", i, " iterations. Difference: ",
                format(difference, digits = 5))
      else
        message("No convergence after ", i, " iterations. Difference: ",
                format(difference, digits = 5))
    }
  }

  if (length(dim(X)) < r + 1)
    X = array(X, dim = c(1, dim(X)))

  if (burnin > 0) {
    if (r == 1)
      inds = list(seq(n) + burnin)
    else
      inds = lapply(lapply(X = n, FUN = seq), "+", burnin)
    X = do.call(`[`, c(list(X), TRUE, inds, drop = FALSE))
  }

  dimnames(X) = list(basisobj$names)

  fsd.fd(X, basisobj)
}
