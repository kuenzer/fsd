#' Analyse functional spatial filters
#'
#' This function is used to analyse functional spatial filters, indicating how
#' the squared norms lie.
#' @param object the filter.
#' @param use.norm which norm to use for the lags.
#' @param make.plot a boolean indicating whether to plot the result.
#' @param ... currently unused
#' @return object list with components \item{balls}{a data frame containing the
#'   cumulative sums of the squared norms of the filter coefficients up to some
#'   maximum lag.} \item{zeroplanes}{a matrix containing the sums of the squared
#'   norms of the filter coefficients lying on the (hyper-) planes with one
#'   dimension of the lag equals zero, along with the sum over all lags as
#'   comparison.} \item{abs.sum}{a data frame containing the cumulative sums of
#'   the norms of the filter coefficients.}
#' @seealso \code{\link{fsd.plot.filters}}
#' @keywords fsd
#' @export
#' @examples
#' \dontrun{
#' summary.fsd.filter(object)
#' }

summary.fsd.filter = function (object, ..., use.norm = Inf, make.plot = FALSE)
{
  npc = dim(object$operators[[1]])[2]
  if (!is.null(object$basis)) {
    B = inprod(object$basis, object$basis)
  } else {
    B = diag(dim(object$operators[[1]])[1])
  }

  if (use.norm < Inf)
    lagnorm = sapply(object$laglist,
                     function(h) {(sum(abs(h)^use.norm))^(1/use.norm)})
  else
    lagnorm = sapply(object$laglist, function(h) {max(abs(h))})

  phinorm = sapply(object$operators, function(M) {diag(t(M) %*% B %*% M)})
  if (npc == 1)
    phinorm = t(phinorm)

  uniquelags = sort(unique(lagnorm))
  normsums = array(dim = c(length(uniquelags), npc),
                   dimnames = list(NULL,
                                   paste0("PC", 1:npc)))

  normsuma = normsums

  for (k in 1:length(uniquelags)) {
    normsums[k, ] = apply(phinorm[ , lagnorm <= uniquelags[k], drop = FALSE], 1,
                          sum)
    normsuma[k, ] = apply(sqrt(phinorm[, lagnorm <= uniquelags[k],
                                       drop = FALSE]), 1, sum)
  }

  normsums2 = sapply(1:length(object$laglist[[1]]),
                     function(i) {apply(phinorm[ , sapply(object$laglist,
                                                          '[', i) == 0,
                                                 drop = FALSE], 1, sum)})
  if (is.vector(normsums2))
    normsums2 = t(normsums2)
  normsums2 = cbind(normsums2, apply(phinorm, 1, sum))

  dimnames(normsums2) = list(paste0("PC", 1:npc),
                             c(paste0("Dim", 1:length(object$laglist[[1]]), 
                                      "=0"), "Total"))

  if (make.plot) {
    matplot(uniquelags, normsums, pch = 1, ylim = 0:1,
            xlab = "Lag", ylab = "Sum of squared norms")
    matplot(uniquelags, normsuma, pch = 1,
            xlab = "Lag", ylab = "Sum of norms")
  }

  list(balls = data.frame(max.lag = uniquelags, normsums),
       zeroplanes = t(normsums2),
       abs.sum = data.frame(max.lag = uniquelags, normsuma))
}
