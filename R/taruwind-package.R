#' TaruWind: reduced-mode terrain and surface wind downscaling
#'
#' New to TaruWind? \href{../doc/getting-started.html}{Click here to open the
#' complete HTML getting-started guide}. It walks through one connected
#' example, from coordinates and downloads to interpreting the result.
#'
#' Start with \code{\link{taruwind_downscale}} for coordinates and a period,
#' or \code{\link{taruwind_build}} and \code{\link{taruwind_predict}} for
#' metric-grid inputs.
#' If the help viewer does not open the link, run
#' \code{browseURL(system.file("doc", "getting-started.html", package = "taruwind"))}.
#' @keywords internal
#' @importFrom Rcpp evalCpp
#' @useDynLib taruwind, .registration = TRUE
"_PACKAGE"
