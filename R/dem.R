#' Download Copernicus GLO-30/GLO-90 DEM tiles from the public AWS registry
#'
#' GLO-90 is the default. If explicitly requested, GLO-30 is attempted first
#' and missing tiles fall back to GLO-90. The product
#' is a digital surface model and its licence/attribution remain the user's
#' responsibility. See <https://registry.opendata.aws/copernicus-dem/>.
#' Requests are limited to 16 source tiles; this helper downloads whole DEM
#' tiles and reuses them from cache. It does not silently substitute ocean
#' zero elevations when a source tile is unavailable.
#'
#' @param bbox Named or unnamed vector `c(xmin, ymin, xmax, ymax)` in degrees.
#' @param target_crs Projected target CRS accepted by `terra`.
#' @param resolution Target resolution in metres.
#' @param cache_dir Existing or creatable directory for source COG tiles.
#' @param prefer_glo30 Whether to try public GLO-30 before universal GLO-90.
#' @param output Optional projected GeoTIFF path; existing files are not replaced.
#' @return A projected `terra::SpatRaster` with source/citation attributes.
#' @export
taruwind_download_copernicus_dem <- function(bbox, target_crs, resolution = 90,
                                        cache_dir = file.path(tempdir(), "taruwind-dem"),
                                        prefer_glo30 = FALSE, output = NULL) {
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("taruwind_download_copernicus_dem requires the optional package 'terra'")
  }
  if (length(bbox) != 4L || any(!is.finite(bbox))) {
    stop("bbox must contain four finite longitude/latitude values")
  }
  names(bbox) <- c("xmin", "ymin", "xmax", "ymax")
  if (bbox["xmin"] >= bbox["xmax"] || bbox["ymin"] >= bbox["ymax"] ||
      bbox["ymin"] < -90 || bbox["ymax"] > 90 ||
      bbox["xmin"] < -180 || bbox["xmax"] > 180) {
    stop("bbox is invalid")
  }
  if(length(resolution)!=1L || !is.finite(resolution) || resolution<=0)
    stop("resolution must be a positive finite scalar")
  if(!is.logical(prefer_glo30)||length(prefer_glo30)!=1L||is.na(prefer_glo30))
    stop("prefer_glo30 must be TRUE or FALSE")
  if(!is.null(output) && file.exists(output)) stop("output already exists")
  tiles_n <- (ceiling(bbox["ymax"])-floor(bbox["ymin"])) *
    (ceiling(bbox["xmax"])-floor(bbox["xmin"]))
  if(tiles_n>16) stop("At most 16 DEM tiles per request; split large requests explicitly")
  if (!dir.exists(cache_dir) && !dir.create(cache_dir, recursive = TRUE)) {
    stop("could not create cache_dir")
  }
  latitude_degrees <- seq.int(floor(bbox["ymin"]), ceiling(bbox["ymax"]) - 1L)
  longitude_degrees <- seq.int(floor(bbox["xmin"]), ceiling(bbox["xmax"]) - 1L)
  files <- character()
  products <- character()
  for (latitude in latitude_degrees) {
    for (longitude in longitude_degrees) {
      tile <- .taruwind_copernicus_tile(latitude, longitude)
      choices <- if (prefer_glo30) c("30m", "90m") else "90m"
      acquired <- FALSE
      for (product in choices) {
        arcseconds <- if (product == "30m") "10" else "30"
        directory <- paste0("Copernicus_DSM_COG_", arcseconds, "_", tile, "_DEM")
        url <- paste0("https://copernicus-dem-", product,
                      ".s3.amazonaws.com/", directory, "/", directory, ".tif")
        destination <- file.path(cache_dir, paste0(directory, ".tif"))
        if (!file.exists(destination)) {
          status <- tryCatch(
            suppressWarnings(utils::download.file(url, destination, mode = "wb", quiet = TRUE)),
            error = function(error) 1L
          )
          if (!identical(status, 0L)) {
            if (file.exists(destination)) unlink(destination)
            next
          }
        }
        files <- c(files, destination)
        products <- c(products, product)
        acquired <- TRUE
        break
      }
      if (!acquired) stop("no public Copernicus DEM tile was found for ", tile)
    }
  }
  merged <- terra::vrt(files)
  cropped <- terra::crop(
    merged,
    terra::ext(bbox["xmin"], bbox["xmax"], bbox["ymin"], bbox["ymax"])
  )
  projected <- terra::project(cropped, target_crs, method = "bilinear", res = resolution)
  if (!is.null(output)) {
    if (file.exists(output)) stop("output already exists")
    terra::writeRaster(projected, output, overwrite = FALSE,
                       gdal = c("COMPRESS=DEFLATE"))
  }
  attr(projected, "taruwind_source_products") <- unique(products)
  attr(projected, "taruwind_citation") <- paste(
    "Copernicus Digital Elevation Model (DEM) was accessed on",
    as.character(Sys.Date()), "from https://registry.opendata.aws/copernicus-dem/."
  )
  projected
}

.taruwind_copernicus_tile <- function(latitude, longitude) {
  latitude_code <- sprintf("%s%02d_00", if (latitude >= 0) "N" else "S", abs(latitude))
  longitude_code <- sprintf("%s%03d_00", if (longitude >= 0) "E" else "W", abs(longitude))
  paste(latitude_code, longitude_code, sep = "_")
}
