#' Creating a functional spatial data object
#'
#' This function extends the possibilities of the function fd found in the
#' package fda to include functional data on a spatial grid of arbitrary
#' dimension r.
#' @param coef an array of arbitrary dimension 1 + r, where the first dimension
#'   corresponds to basis functions.
#' @param basisobj a functional basis object defining the basis.
#' @param fdnames a list of length r with each member containing the labels for
#'   the dimensions of the data.
#' @return the \code{fsd.fd} object.
#' @keywords fsd
#' @import fda
#' @export
#' @examples
#' \dontrun{
#' fsd.fd(coef, basisobj)
#' }

fsd.fd = function (coef = NULL, basisobj = NULL, fdnames = NULL)
{
  if (is.null(coef) && is.null(basisobj))
    basisobj <- create.bspline.basis(nbasis = 2, norder = 1)
  if (is.null(coef))
    coef <- rep(0, basisobj[["nbasis"]])
  type <- basisobj$type
  {
    if (!is.numeric(coef))
      stop("'coef' is not numeric.")
    else if (is.vector(coef)) {
      coef <- as.matrix(coef)
      if (identical(type, "constant"))
        coef <- t(coef)
      coefd <- dim(coef)
      ndim <- length(coefd)
    }
    else if (is.array(coef)) {
      coefd <- dim(coef)
      ndim <- length(coefd)
    }
    else stop("Type of 'coef' is not correct")
  }
  if (is.null(basisobj)) {
    dimC <- dim(coef)
    nb <- {
      if (is.null(dimC))
        length(coef)
      else dimC[1]
    }
    basisobj <- create.fourier.basis(nbasis = nb)
    type <- basisobj$type
  }
  else if (!(inherits(basisobj, "basisfd")))
    stop("Argument basis must be of basis class")
  nbasis <- basisobj$nbasis
  dropind <- basisobj$dropind
  ndropind <- length(basisobj$dropind)
  if (coefd[1] != nbasis - ndropind)
    stop("First dim. of 'coef' not equal to 'nbasis - ndropind'.")
  if (ndim > 1)
    nrep <- coefd[2]
  else nrep <- 1
  if (ndim > 2)
    nvar <- coefd[3]
  else nvar <- 1
  if (is.null(fdnames)) {
    fdnames = NULL
    if (!is.vector(coef)) {
      fdnames = dimnames(coef)
    }
    if (is.null(fdnames)) {
      fdnames <- c("time",
                   lapply(1:(ndim-1),
                          function(k) {paste0("Dim", k, " ",
                                              as.character(1:coefd[1+k]))}))
    }
    for (i in 2:length(fdnames) - 1) {
      if (is.null(fdnames[[i+1]]))
        fdnames[[i+1]] = paste0("Dim", i, " ", as.character(1:coefd[1+i]))
    }
    if (is.null(names(fdnames))) {
      names(fdnames) <- c("args", paste0("Dim", as.character(1:(ndim-1))))
    }
  }
  if (is.null(dimnames(coef))) {
    dimc <- dim(coef)
    ndim <- length(dimc)
    dnms <- vector("list", ndim)
    rdms <- mapply(function(a,b) {a == length(b)}, dimc, fdnames)
    dnms[[1]] = basisobj$names
    for (i in which(rdms)) {
      dnms[[i]] <- fdnames[[i]]
    }
    dimnames(coef) <- dnms
  }
  fdobj <- list(coefs = coef, basis = basisobj, fdnames = fdnames)
  oldClass(fdobj) <- c("fsd.fd", "fd")
  fdobj
}


#' Converting a functional data object into a functional spatial data object
#'
#' Works as one would expect
#' @param x a functional data object or something that can be converted into
#'   one.
#' @param ... currently unused
#' @seealso \code{\link{fsd.fd}}
#' @keywords fsd
#' @export
#' @examples
#' \dontrun{
#' as.fsd.fd(x)
#' }

as.fsd.fd = function (x, ...)
{
  if (!inherits(x, "fd"))
    x = as.fd(x)

  r = length(dim(x$coefs)) - 1

  fdnames = x$fdnames[1:(r+1)]

  fsd.fd(coef = x$coefs, basisobj = x$basis, fdnames = fdnames)
}


#' Converting a functional spatial data object into a functional data object
#'
#' Works as one would expect
#' @param x a functional spatial data object
#' @param ... currently unused
#' @seealso \code{\link{fsd.fd}}
#' @keywords fsd
#' @import fda
#' @export
#' @examples
#' \dontrun{
#' as.fd(x)
#' }

as.fd.fsd.fd = function (x, ...)
{
  if (!inherits(x, "fsd.fd"))
    stop("x needs to be a functional spatial data object")

  x = x[drop = TRUE]

  r = length(dim(x$coefs)) - 1

  if (r > 2)
    stop("x can only have a two-dimensional spatial grid")

  if (length(x$fdnames) == 2)
    x$fdnames = c(x$fdnames, list(funs = "values"))
  if (length(x$fdnames) > 3)
    x$fdnames = x$fdnames[1:3]

  fd(x$coefs, basisobj = x$basis, fdnames = x$fdnames)
}


#' Subsetting a functional spatial data object
#'
#' Subsetting works as with conventional arrays. By default, dimensions can be 
#' dropped.
#' @param fsdobj a functional spatial data object
#' @param ... the indices to take
#' @param drop whether to flatten the spatial grid by dropping dimensions of
#'   length one.
#' @keywords fsd
#' @export
#' @examples
#' data("temp")
#' plot(temp[1:2, 1:2, 1, drop = TRUE])

`[.fsd.fd` = function (fsdobj, ..., drop = TRUE)
{
  r = length(dim(fsdobj$coefs)) - 1

  args = lapply( eval(substitute(alist(...))) ,
                 function(y) {if (is.call(y))
                                eval.parent(y, n = 3)
                              else {
                                if (length(y) > 1)
                                  y
                                else if (nchar(as.character(y)) == 0)
                                  y
                                else
                                  eval.parent(y, n = 3)
                              }})

  inds = c(TRUE, args)
  if (length(inds) == 1)
    inds = c(inds, rep(TRUE, r))

  fsdobj$coefs = do.call(what = '[', args = c(list(fsdobj$coefs),
                                              inds, drop = drop))

  fsdobj$fdnames = fsdobj$fdnames[1:(r+1)]

  fsdobj$fdnames = mapply('[', fsdobj$fdnames, inds, SIMPLIFY = FALSE)

  if (drop) {
    collapsing = (sapply(fsdobj$fdnames, length) == 1)
    collapsing[1] = FALSE
    if (is.vector(fsdobj$coefs)) {
      collapsing[2] = FALSE
      fsdobj$coefs = as.matrix(fsdobj$coefs)
    }
    fsdobj$fdnames = fsdobj$fdnames[!collapsing]
  }

  fsdobj
}


#' Subsetting a functional spatial data object
#'
#' Take only every \code{s}-th element
#' @param fsdobj a functional spatial data object
#' @param s an integer or a vector with one entry for every dimension of the
#'   spatial grid.
#' @keywords fsd
#' @export
#' @examples
#' data("temp")
#' plot(temp %/% 3)

`%/%.fsd.fd` = function (fsdobj, s)
{
  if(!inherits(fsdobj, "fsd.fd"))
    stop("Needs functional spatial data")
  if (length(s) > 1 && length(s) != length(dim(fsdobj$coefs)) -1 )
    stop("Incompatible dimensions")
  if (any(s < 1) || any(s %% 1 != 0))
    stop("Needs positive integer s")

  r = length(dim(fsdobj$coefs)) - 1

  if (length(s) == 1)
    s = rep(s, r)

  inds = lapply(1:r, function(i) {seq(from = 1, to = dim(fsdobj$coefs)[1+i],
                                      by = s[i])})

  do.call('[', c(list(fsdobj), inds, drop = FALSE))
}


#' Arithmetic for functional spatial data objects
#'
#' Works as one would expect.
#' @param fsd1,fsd2 functional spatial data objects
#' @keywords fsd
#' @export
#' @examples
#' \dontrun{
#' fsdobj1 + fsdobj2
#' }

`+.fsd.fd` = function (fsd1, fsd2)
{
  if (!inherits(fsd1, "fsd.fd") || !inherits(fsd2, "fsd.fd"))
    stop("Needs compatible data")
  else if (!(fsd1$basis == fsd2$basis))
    stop("Incompatible bases")

  if (length(dim(fsd1$coefs)) != length(dim(fsd2$coefs))) {
    if (length(dim(fsd2$coefs)) == 2 && dim(fsd2$coefs)[2] == 1)
      fsd1$coefs = fsd1$coefs + as.vector(fsd2$coefs)
    else if (length(dim(fsd1$coefs)) == 2 && dim(fsd1$coefs)[2] == 1) {
      tmp = as.vector(fsd1$coefs)
      fsd1 = fsd2
      fsd1$coefs = fsd1$coefs + tmp
    } else
      stop("Incompatible dimensions")
  }
  else if (any(dim(fsd1$coefs) != dim(fsd2$coefs)))
    stop("Incompatible dimensions")
  else
    fsd1$coefs = fsd1$coefs + fsd2$coefs

  fsd1
}


#' Arithmetic for functional spatial data objects
#'
#' Works as one would expect.
#' @param fsd1,fsd2 functional spatial data objects
#' @keywords fsd
#' @export
#' @examples
#' \dontrun{
#' fsdobj1 - fsdobj2
#' - fsdobj
#' }

`-.fsd.fd` = function (fsd1, fsd2 = NULL)
{
  if (is.null(fsd2)) {
    fsd1$coefs = -fsd1$coefs
    fsd1
  } else {
    fsd2 = -fsd2
    fsd1 + fsd2
  }
}


#' Arithmetic for functional spatial data objects
#'
#' Works as one would expect.
#' @param scalar a scalar number
#' @param fsdobj a functional spatial data object
#' @keywords fsd
#' @export
#' @examples
#' \dontrun{
#' 2 * fsdobj
#' }

`*.fsd.fd` = function (scalar, fsdobj)
{
  if (is.numeric(fsdobj)) {
    tmp = fsdobj
    fsdobj = scalar
    scalar = tmp
  }

  if (!is.numeric(scalar) ||
      length(scalar) != 1 ||
      !inherits(fsdobj, "fsd.fd"))
    stop("Only scalar multiplication with an fsd.fd object is possible")

  fsdobj$coefs = scalar * fsdobj$coefs

  fsdobj
}


#' Arithmetic for functional spatial data objects
#'
#' Works as one would expect.
#' @param fsdobj a functional spatial data object
#' @param scalar a scalar number
#' @keywords fsd
#' @export
#' @examples
#' \dontrun{
#' fsdobj / 2
#' }

`/.fsd.fd` = function (fsdobj, scalar)
{
  if (!is.numeric(scalar) ||
      length(scalar) != 1 ||
      !inherits(fsdobj, "fsd.fd"))
    stop("Only scalar division of an fsd.fd object is possible")

  fsdobj$coefs = fsdobj$coefs / scalar

  fsdobj
}

