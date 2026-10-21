#' Automatic-workflow safety limits
#'
#' TaruWind's automatic download workflow is intended for point networks and
#' local terrain domains, not accidental continental raster production. These
#' conservative limits are checked before network or DEM work starts. Advanced
#' users of [taruwind_download_era5()] can explicitly opt out with
#' `allow_large = TRUE`. The main [taruwind_downscale()] workflow does not expose
#' this override; split larger requests into recoverable batches.
#' A point or raw ERA5 request is limited to ten calendar years per call.
#' Local area maps are limited to 168 hourly timestamps (seven days). Longer
#' series can be built in batches with reusable prepared responses. Terrain-cell
#' limits count the reflected computational grid, not just the output square.
#' They do not estimate network bytes or guarantee that a large request will fit
#' in memory.
#'
#' @return A named list containing the default limits.
#' @export
taruwind_limits <- function() {
  list(
    supported_era5_wind_heights_m = c(10L, 100L),
    max_sites = 100L,
    max_site_hours = 1e7,
    max_years_per_request = 10L,
    max_area_hours = 168L,
    max_cells_per_site = 2e6,
    max_total_terrain_cells = 6e7
  )
}

.taruwind_request_hours <- function(period) {
  as.numeric(difftime(period$end, period$start, units = "hours")) + 1
}

.taruwind_calendar_year_limit <- function(start, years) {
  boundary <- as.POSIXlt(start, tz = "UTC")
  boundary$year <- boundary$year + years
  as.POSIXct(boundary, tz = "UTC")
}

.taruwind_guard_request <- function(
    sites, period, resolution = NULL, domain_km = NULL,
    allow_large = FALSE, allow_override = TRUE) {
  if (!is.logical(allow_large) || length(allow_large) != 1L || is.na(allow_large)) {
    stop("allow_large must be TRUE or FALSE")
  }
  limits <- taruwind_limits()
  hours <- .taruwind_request_hours(period)
  site_hours <- nrow(sites) * hours
  issues <- character()
  if (nrow(sites) > limits$max_sites) {
    issues <- c(issues, paste0(nrow(sites), " sites > ", limits$max_sites))
  }
  if (site_hours > limits$max_site_hours) {
    issues <- c(
      issues,
      paste0(format(site_hours, scientific = FALSE, trim = TRUE),
             " site-hours > ",
             format(limits$max_site_hours, scientific = FALSE, trim = TRUE))
    )
  }
  if (period$end >= .taruwind_calendar_year_limit(
      period$start, limits$max_years_per_request)) {
    issues <- c(issues, paste0("requested period exceeds ",
      limits$max_years_per_request, " calendar years per call"))
  }
  if (!is.null(resolution) && !is.null(domain_km)) {
    if (!is.numeric(resolution) || length(resolution) != 1L ||
        !is.finite(resolution) || resolution <= 0 ||
        !is.numeric(domain_km) || length(domain_km) != 1L ||
        !is.finite(domain_km) || domain_km <= 0) {
      stop("resolution and domain_km must be positive finite scalars", call. = FALSE)
    }
    n <- max(3, 2 * ceiling((domain_km * 1000 / resolution - 1) / 2) + 1)
    cells_per_site <- 9 * n * n
    total_cells <- nrow(sites) * cells_per_site
    if (cells_per_site > limits$max_cells_per_site) {
      issues <- c(
        issues,
        paste0(format(round(cells_per_site), big.mark = ",", scientific = FALSE),
               " reflected computational cells/site > ",
               format(limits$max_cells_per_site, big.mark = ",", scientific = FALSE))
      )
    }
    if (total_cells > limits$max_total_terrain_cells) {
      issues <- c(
        issues,
        paste0(format(round(total_cells), big.mark = ",", scientific = FALSE),
               " total reflected computational cells > ",
               format(limits$max_total_terrain_cells, big.mark = ",", scientific = FALSE))
      )
    }
  }
  if (!length(issues)) return(invisible(TRUE))
  message <- paste0(
    "request exceeds TaruWind automatic-workflow safety limits: ",
    paste(issues, collapse = "; "),
    if (allow_override) {
      ". Reduce/split the request or set allow_large = TRUE deliberately."
    } else {
      ". Reduce or split the request into recoverable batches."
    }
  )
  if (!allow_large) stop(message, call. = FALSE)
  warning(message, call. = FALSE)
  invisible(TRUE)
}
