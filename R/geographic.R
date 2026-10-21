#' Prepare a geographic target using raster terrain and roughness
#'
#' A local azimuthal-equidistant grid is centred on the target, so its axes at
#' the target coincide with east/north wind components. Raster rows are reversed
#' explicitly for the native solver. No missing cells are filled implicitly.
#' @param lon,lat One target longitude and latitude in degrees.
#' @param roughness Aerodynamic roughness SpatRaster or raster filename, in
#'   metres, or the explicit string "worldcover" for automatic acquisition.
#'   Land-cover codes must otherwise be converted before calling this function.
#' @param depth,source_z0 Explicit neutral reference depth and effective source
#'   roughness in metres; these are not inferred from target land cover.
#' @param dem Optional elevation SpatRaster or filename. If NULL, GLO-90 is
#'   downloaded with [taruwind_download_copernicus_dem()].
#' @param resolution Metric calculation spacing, default 90 metres.
#' @param domain_m Requested square domain width in metres, rounded upwards to
#'   an odd number of cells with the target at the central cell centre.
#' @param height Evaluation height, 10 or 100 metres.
#' @param cache_dir Directory for downloaded DEM tiles.
#' @param quiet Suppress preparation messages.
#' @return A prepared response usable with [taruwind_predict()], with geographic
#'   input provenance. Caller-provided rasters must have a valid CRS.
#' @export
taruwind_prepare_site <- function(lon, lat, roughness, depth, source_z0,
    dem=NULL, resolution=90, domain_m=30000, height=10,
    cache_dir=file.path(file.path(tempdir(), "taruwind-cache"), "dem"), quiet=TRUE) {
  q <- .taruwind_site_inputs(lon,lat,roughness,depth,source_z0,dem,
    resolution,domain_m,height,cache_dir,quiet)
  op <- do.call(taruwind_build,q$arguments)
  op$provenance$geographic <- q$geographic
  op
}

.taruwind_site_inputs <- function(lon, lat, roughness, depth, source_z0,
    dem=NULL, resolution=90, domain_m=30000, height=10,
    cache_dir=file.path(file.path(tempdir(), "taruwind-cache"), "dem"), quiet=TRUE,
    return_layers=FALSE) {
  .taruwind_require("terra", "geographic preparation")
  target <- .taruwind_sites(lon=lon,lat=lat)
  if(nrow(target)!=1L || abs(lat)>85) stop("One target between 85 S and 85 N required")
  height <- .taruwind_wind_height(height)
  if(length(resolution)!=1L || !is.finite(resolution) || resolution<=0 ||
     length(domain_m)!=1L || !is.finite(domain_m) || domain_m<=0)
    stop("Positive finite resolution and domain_m required")
  n <- max(3, 2*ceiling((domain_m/resolution-1)/2)+1)
  if(9*n*n>2e6) stop("Requested domain exceeds the reflected-cell limit")
  crs <- sprintf("+proj=aeqd +lat_0=%.10f +lon_0=%.10f +datum=WGS84 +units=m",lat,lon)
  half <- n*resolution/2
  grid <- terra::rast(nrows=n,ncols=n,xmin=-half,xmax=half,
                       ymin=-half,ymax=half,crs=crs)
  read_layer <- function(x) {
    r <- if(inherits(x,"SpatRaster")) x else terra::rast(x)
    if(terra::nlyr(r)!=1L || !nzchar(terra::crs(r))) stop("Single-layer raster with CRS required")
    r
  }
  if(is.null(dem)) {
    if (!quiet) message("  Terrain: acquiring GLO-90 DEM (cached tiles are reused).")
    # Densified boundary, including a halo for reprojection interpolation.
    a <- seq(-half-resolution*2,half+resolution*2,length.out=33)
    b <- half+resolution*2
    edge <- rbind(cbind(a,-b),cbind(a,b),cbind(-b,a),cbind(b,a))
    ll <- terra::crds(terra::project(terra::vect(edge,crs=crs),"EPSG:4326"))
    bbox <- c(min(ll[,1]),min(ll[,2]),max(ll[,1]),max(ll[,2]))
    if(diff(range(ll[,1]))>180) stop("Antimeridian-crossing DEM downloads are not supported")
    dem <- taruwind_download_copernicus_dem(bbox,crs,resolution,cache_dir)
  }
  align <- function(r,method) {
    if(terra::compareGeom(r,grid,stopOnError=FALSE)) r else
      terra::project(r,grid,method=method)
  }
  h <- align(read_layer(dem),"bilinear")
  # Nearest preserves supplied aerodynamic classes instead of inventing a
  # physically unsupported arithmetic interpolation between roughness lengths.
  automatic_roughness <- identical(roughness,"worldcover")
  if (!quiet) message(if (automatic_roughness)
    "  Roughness: reading WorldCover windows or cached raster windows." else
    "  Roughness: aligning caller-supplied aerodynamic roughness.")
  z <- if(automatic_roughness) taruwind_download_worldcover(grid,
    cache_dir=file.path(cache_dir,"worldcover")) else
    align(read_layer(roughness),"near")
  south_first <- function(r) {
    m <- terra::as.matrix(r,wide=TRUE)
    m[nrow(m):1L,,drop=FALSE]
  }
  arguments <- list(dem=south_first(h),roughness=south_first(z),dx=resolution,
    x0=-half+resolution/2,y0=-half+resolution/2,x=0,y=0,depth=depth,
    source_z0=source_z0,height=height)
  geographic <- list(lon=lon,lat=lat,crs=crs,resolution=resolution,
    actual_domain_m=n*resolution,dem_resampling="bilinear",roughness_resampling="nearest",
    wind_frame="east/north at target",roughness_source=if(automatic_roughness)
      attr(z,"taruwind_roughness_provenance") else "caller supplied")
  if(return_layers) {
    names(h) <- "elevation_m"
    names(z) <- "roughness_m"
  }
  list(arguments=arguments,geographic=geographic,
    layers=if(return_layers) c(h,z) else NULL)
}

#' Retrieve the terrain layers used for one geographic site
#'
#' Returns the GLO-90 elevation and aerodynamic roughness after alignment to
#' the same target-centred metric grid used by [taruwind_prepare_site()]. The
#' layers can be plotted with [taruwind_plot_map()]. No wind is downloaded or
#' response constructed. WorldCover roughness is a land-cover lookup, not a
#' measured roughness map.
#' @inheritParams taruwind_prepare_site
#' @return A two-layer `terra::SpatRaster` named `elevation_m` and `roughness_m`.
#' @export
taruwind_site_layers <- function(lon,lat,roughness="worldcover",dem=NULL,
    resolution=90,domain_m=30000,
    cache_dir=file.path(file.path(tempdir(),"taruwind-cache"),"dem"),quiet=TRUE) {
  .taruwind_site_inputs(lon,lat,roughness,depth=NULL,source_z0=NULL,
    dem=dem,resolution=resolution,domain_m=domain_m,height=10,
    cache_dir=cache_dir,quiet=quiet,return_layers=TRUE)$layers
}

#' Download ERA5 and downscale geographic target time series
#'
#' Downloads ERA5 and (unless supplied) GLO-90. Roughness and source-profile
#' parameters remain explicit: this function does not equate source roughness
#' with target roughness or assume that ERA5 has one global roughness value.
#' Downloads and raster preparation are serial to avoid cache races. Independent
#' numeric response construction uses [taruwind_build_many()] with its memory
#' guard and portable workers. Prepared responses can then be reused.
#' This is a point workflow, not a national atlas tool.
#' @inheritParams taruwind_download_era5
#' @param roughness Raster or filename shared by all sites, or a list of one
#'   raster/filename per site, or "worldcover" for automatic acquisition.
#'   Raster values must be aerodynamic roughness in metres.
#' @param depth,source_z0 Scalar or one value per site, in metres. If source_z0
#'   is NULL, [taruwind_source_profile()] retrieves the 2021 reference FSR and
#'   applies the documented diagnostic convention. If depth is NULL it is
#'   computed by [taruwind_neutral_depth()]. Explicit inputs override automation.
#' @param dem NULL for automatic GLO-90, one raster/filename, or a per-site list.
#' @param resolution,domain_m Grid spacing and local domain width in metres.
#' @param responses Optional named list of prepared geographic responses from
#'   a previous call, indexed by site_id. Coordinates and heights are checked.
#' @param workers "auto" (half the detected cores, capped at 16 and site count),
#'   or an integer from 1 to 16; the memory guard may reduce concurrency.
#' @param profile_year Reference FSR climatology year; default 2021, not inferred
#'   from the wind period. See [taruwind_source_profile()].
#' @param profile_wait_seconds CDS queue polling limit per site/call.
#' @return A tibble (also a data frame) with one row per site-hour.
#'   `wind_speed` is vector-bilinear ERA5 speed; `taruwind_speed` is the
#'   downscaled speed. `attr(result, "responses")` is a named list of prepared
#'   responses, one per site. `attr(result, "forcing_provenance")` records ERA5
#'   source and interpolation metadata. These two attributes are separate
#'   from the printed columns. The current ERA5 download metadata is also
#'   forwarded as attributes, including `sites`, `source`, `license`,
#'   `citation`, `retrieved_utc`, `wind_height_m`, `component_variables`,
#'   `interpolation` and `cache_dir`. Inspect `names(attributes(result))` for
#'   the complete set. `saveRDS()` retains attributes; CSV does not.
#' @export
taruwind_downscale <- function(sites=NULL,lon=NULL,lat=NULL,start=NULL,end=NULL,
    period=NULL,roughness=NULL,depth=NULL,source_z0=NULL,dem=NULL,
    resolution=90,domain_m=30000,wind_height=10,local_tz=FALSE,
    cache_dir=file.path(tempdir(), "taruwind-cache"),api_key=NULL,
    responses=NULL,workers="auto",profile_year=2021,profile_wait_seconds=60,
    refresh=c("auto","always","never"),quiet=FALSE) {
  if (!is.logical(quiet) || length(quiet)!=1L || is.na(quiet)) stop("quiet must be TRUE or FALSE")
  refresh <- match.arg(refresh)
  report = function(...) if (!quiet) message(...)
  target <- .taruwind_sites(sites,lon,lat)
  when <- .taruwind_resolve_period(start,end,period)
  if (!is.numeric(resolution) || length(resolution)!=1L ||
      !is.finite(resolution) || resolution<=0 ||
      !is.numeric(domain_m) || length(domain_m)!=1L ||
      !is.finite(domain_m) || domain_m<=0) {
    stop("resolution and domain_m must be positive finite scalars", call.=FALSE)
  }
  .taruwind_guard_request(target,when,
    resolution=if(is.null(responses)) resolution else NULL,
    domain_km=if(is.null(responses)) domain_m/1000 else NULL,
    allow_override=FALSE)
  wind_height <- .taruwind_wind_height(wind_height)
  if(is.null(responses) && is.null(roughness)) {
    stop("Supply roughness or explicitly select roughness = 'worldcover'")
  }
  hours <- .taruwind_request_hours(when)
  report("Request plan: ", nrow(target), " site(s) x ",
    format(hours,big.mark=",",scientific=FALSE), " hours = ",
    format(nrow(target)*hours,big.mark=",",scientific=FALSE),
    " site-hours. Network bytes depend on source chunking and cache state.")
  report("Cache directory: ", normalizePath(cache_dir, mustWork=FALSE))
  report("[1/6] Checking currently available ERA5 hours before terrain work.")
  axis <- .taruwind_era5_axis(when,api_key,file.path(cache_dir,"era5"),refresh,quiet)
  report("[2/6] Source profile: ", if (!is.null(responses)) "reusing prepared responses" else
    if (!is.null(source_z0)) "using supplied values" else
    "retrieving reference FSR (cache first; CDS queue may require waiting)")
  if(is.null(responses)) {
    profiles <- NULL
    if(is.null(source_z0)) {
      profiles <- lapply(seq_len(nrow(target)),function(i)
        taruwind_source_profile(target$lon[i],target$lat[i],wind_height,
          reference_year=profile_year,api_key=api_key,
          cache_dir=file.path(cache_dir,"source_profile"),wait_seconds=profile_wait_seconds,quiet=quiet))
      source_z0 <- vapply(profiles,`[[`,0.,"source_z0")
    }
    if(is.null(depth)) depth <- taruwind_neutral_depth(target$lat)
    per_site <- function(x,i) {
      if(is.list(x)) {
        if(length(x)!=nrow(target)) stop("Per-site list must match target count")
        return(x[[i]])
      }
      x
    }
    scalar <- function(x,i) {
      if(!is.numeric(x)||!length(x)%in%c(1L,nrow(target))) stop("Profile input length must be one or site count")
      x[if(length(x)==1L) 1L else i]
    }
    report("[3/6] Preparing terrain and roughness for ", nrow(target), " site(s).")
    inputs <- lapply(seq_len(nrow(target)),function(i) {
      report("  Site ", i, "/", nrow(target), ": ", target$site_id[i])
      .taruwind_site_inputs(target$lon[i],target$lat[i],per_site(roughness,i),
        scalar(depth,i),scalar(source_z0,i),dem=per_site(dem,i),
        resolution=resolution,domain_m=domain_m,height=wind_height,
        cache_dir=file.path(cache_dir,"dem"),quiet=quiet)
    })
    report("[4/6] Building terrain responses. This step is independent of requested hours.")
    responses <- taruwind_build_many(lapply(inputs,`[[`,"arguments"),workers=workers)
    for(i in seq_along(responses)) responses[[i]]$provenance$geographic <- inputs[[i]]$geographic
    if(!is.null(profiles)) for(i in seq_along(responses))
      responses[[i]]$provenance$source_profile <- profiles[[i]]$provenance
    names(responses) <- target$site_id
  } else {
    report("[3/6] Terrain preparation skipped: reusing responses.")
    report("[4/6] Terrain solve skipped: reusing responses.")
  }
  for(i in seq_len(nrow(target))) {
    op <- responses[[target$site_id[i]]]
    g <- op$provenance$geographic
    if(!inherits(op,"taruwind_modal_response")||is.null(g)||
       abs(g$lon-target$lon[i])>1e-8||abs(g$lat-target$lat[i])>1e-8||
       op$provenance$height!=wind_height) stop("Prepared response coordinates/height do not match")
  }
  report("[5/6] Reading hourly ERA5 wind from the CDS Zarr store.")
  winds <- .taruwind_download_era5_ready(target,when,axis,local_tz,FALSE,
    quiet,wind_height)
  report("[6/6] Applying responses to ", nrow(winds), " site-hour records.")
  out <- tibble::as_tibble(winds)
  out$taruwind_speed <- NA_real_
  for(id in target$site_id) {
    at <- which(out$site_id==id)
    out$taruwind_speed[at] <- taruwind_predict(responses[[id]],out$u[at],out$v[at])$speed
  }
  attr(out,"responses") <- responses
  attr(out,"forcing_provenance") <- attributes(winds)[c("source","license","citation","retrieved_utc","interpolation")]
  report("Completed: ", nrow(out), " rows. Cache retained for reuse.")
  out
}
