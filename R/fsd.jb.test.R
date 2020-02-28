#' Perform a Jarque-Bera Test on Normality of Functional Spatial Data
#'
#' This function performs a test on normality of some functional spatial data.
#' @param X.spca a spectral principal components analysis as performed by
#'   \code{\link{fsd.spca}}, i.e. a list that in particular contains the
#'   estimated SPC scores and the estimated spectral density.
#' @param Npc the number of spectral principal components to use for the test.
#' @param L the maximum lag to compute the covariance of the scores during the
#'   computation of the autocovariances.
#' @param var.method the method that is used to calculate the long-run variance
#'   of the SFPC scores. Either "direct" or "integral".
#' @return A list with components \item{T4}{the test statistic.}
#'   \item{p.value}{the p-value.} \item{df}{the degrees of freedom.}
#'   \item{T4.vector}{the vector of the test statistics for each single SPC.}
#' @seealso \code{\link{fsd.spca}}, \code{\link{fsd.spca.cov}}
#' @details To ensure accuracy of numerical integration during the Fourier
#'   transform, the frequencies of \code{F} should be a suitably dense grid.
#' @keywords fsd
#' @import stats
#' @export
#' @examples
#' \dontrun{
#' fsd.jb.test(F, L)
#' }

fsd.jb.test = function (X.spca, Npc = NULL, L = NULL, var.method = "integral")
{
  if (is.null(X.spca$F) && var.method != "direct")
    stop("X.spca needs to contain the spectral density.")
  if (is.null(X.spca$scores))
    stop("X.spca needs to contain the SPC scores.")

  Y = X.spca$scores
  dims = dim(Y)[-1]

  if (is.null(Npc))
    Npc = dim(Y)[1]
  else
    Npc = min(Npc, dim(Y)[1])

  Y = Y - apply(Y, 1, mean)

  mu2 = apply(Y^2, 1, mean)[1:Npc]
  mu3 = apply(Y^3, 1, mean)[1:Npc]
  mu4 = apply(Y^4, 1, mean)[1:Npc]

  if (is.null(L))
    L = dims - 1
  else if (length(L) == 1)
    L = rep(L,length(dims))
  if (var.method == "direct") {
    hlist = fsd:::unfold(lapply(L, function(h) {-h:h}))
    Clist = sapply(hlist, FUN = fsd.covariance, Y = Y, centered = TRUE,
                   simplify = "array")
    if (is.null(dim(Clist))) {
      C = Clist
      dim(C) = c(1, length(C))
    } else
      C = apply(Clist, 3, diag)
  } else
    C = fsd.spca.cov(F = X.spca$F, L = L)$cov

  F3 = apply(C[1:Npc, , drop = FALSE]^3, 1, sum)
  F4 = apply(C[1:Npc, , drop = FALSE]^4, 1, sum)

  T4v = prod(dims) * ( mu3^2 / (6 *F3) + (mu4 - 3*mu2^2)^2 / (24 *F4)  )

  T4 = sum(T4v)

  list(T4 = T4, p.value = pchisq(T4, df = 2*length(T4v), lower.tail = FALSE),
       df = 2*length(T4v), T4.vector = T4v)
}
