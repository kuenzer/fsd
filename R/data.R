#' Temperature Data in Wyoming
#'
#' CPC daily Temperature data of mean temperatures in Wyoming
#' projected onto a B-spline basis of 36 functions for each year,
#' from 1979 to 2017 and with a spatial resolution of 0.5 degrees.
#'
#' The gridded data has been obtained by Shepard's method.
#'
#' @docType data
#'
#' @usage data(temp)
#'
#' @format An object of class \code{"fsd.fd"}.
#'
#' @keywords datasets
#'
#' @references ESRL Physical Sciences Division
#' (\url{https://www.esrl.noaa.gov/psd/})
#'
#' @source \url{https://www.esrl.noaa.gov/psd/data/gridded/data.cpc.globaltemp.html}
#'
#' @examples
#' data(temp)
#' plot(temp)
"temp"
