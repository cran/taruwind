#' Downscale a small local grid of target points
#'
#' Each output cell centre receives its own complete TaruWind response. This
#' produces a map of point estimates on a regular metric grid; it is not a
#' single shared-domain flow solve. Work therefore grows with cell count.
#' The square output grid is centred on `lon, lat` and rounded to an odd
#' number of cells so that one cell represents the centre exactly. Its spacing
#' (`map_resolution`) is independent of the terrain grid spacing (`resolution`).
#' A call may cover at most 168 hourly timestamps (seven complete UTC days).
#' For longer periods, make consecutive calls and reuse the prepared responses
#' from `attr(previous$wind, "responses")`.
#'
#' @param lon,lat Centre coordinates in decimal degrees.
#' @param start,end,period UTC period, as in [taruwind_downscale()].
#' @param area_m Requested output-map width in metres.
#' @param map_resolution Output-cell spacing in metres; does not change the
#'   terrain solver's grid resolution.
#' @param max_cells Maximum output cells, 1 to 100; defaults to 25. Raise this
#'   only when the cost of one response per cell is acceptable.
#' @inheritParams taruwind_downscale
#' @return A list with `wind` (site-hour tibble), `maps` (time-labelled
#'   `terra::SpatRaster` layers for TaruWind and vector-bilinear ERA5), and
#'   `grid` (output grid and coordinate provenance). Use [taruwind_plot_map()]
#'   to display an individual map layer.
#' @export
taruwind_downscale_area <- function(lon,lat,start=NULL,end=NULL,period=NULL,
    area_m=990,map_resolution=500,max_cells=25,
    roughness="worldcover",depth=NULL,source_z0=NULL,dem=NULL,
    resolution=90,domain_m=30000,wind_height=10,
    cache_dir=file.path(tempdir(),"taruwind-cache"),api_key=NULL,
    workers="auto",responses=NULL,profile_year=2021,profile_wait_seconds=60,
    refresh=c("auto","always","never"),quiet=FALSE) {
  .taruwind_require("terra","area downscaling")
  centre <- .taruwind_sites(lon=lon,lat=lat)
  if(nrow(centre)!=1L || abs(lat)>85) stop("One centre between 85 S and 85 N required")
  if(!is.numeric(area_m)||length(area_m)!=1L||!is.finite(area_m)||area_m<=0||
     !is.numeric(map_resolution)||length(map_resolution)!=1L||
     !is.finite(map_resolution)||map_resolution<=0)
    stop("area_m and map_resolution must be positive finite scalars")
  if(!is.numeric(max_cells)||length(max_cells)!=1L||!is.finite(max_cells)||
     max_cells<1||max_cells>100||max_cells!=floor(max_cells))
    stop("max_cells must be an integer between 1 and 100")
  n <- max(1,2*ceiling((area_m/map_resolution-1)/2)+1)
  if(n*n>max_cells) stop("Output map has ",n*n,
    " cells; increase map_resolution, reduce area_m, or explicitly raise max_cells")
  when <- .taruwind_resolve_period(start,end,period)
  hours <- .taruwind_request_hours(when)
  max_hours <- taruwind_limits()$max_area_hours
  if(hours>max_hours) stop("Local area maps are limited to ",max_hours,
    " hourly timestamps (seven days) per call; requested ",hours,
    ". Split the period and reuse prepared responses from a previous area result.",
    call.=FALSE)
  local_crs <- sprintf("+proj=aeqd +lat_0=%.10f +lon_0=%.10f +datum=WGS84 +units=m",lat,lon)
  half <- n*map_resolution/2
  template <- terra::rast(nrows=n,ncols=n,xmin=-half,xmax=half,
    ymin=-half,ymax=half,crs=local_crs)
  xy <- terra::xyFromCell(template,seq_len(terra::ncell(template)))
  coordinates <- terra::crds(terra::project(terra::vect(xy,crs=local_crs),"EPSG:4326"))
  sites <- data.frame(site_id=sprintf("cell_%03d",seq_len(nrow(coordinates))),
    lon=coordinates[,1],lat=coordinates[,2])
  wind <- taruwind_downscale(sites=sites,start=when$start,end=when$end,
    roughness=roughness,depth=depth,source_z0=source_z0,dem=dem,
    resolution=resolution,domain_m=domain_m,wind_height=wind_height,
    cache_dir=cache_dir,api_key=api_key,workers=workers,responses=responses,
    profile_year=profile_year,profile_wait_seconds=profile_wait_seconds,
    refresh=refresh,quiet=quiet)
  timestamps <- sort(unique(wind$timestamp_utc))
  variables <- c("taruwind_speed","wind_speed")
  prefixes <- c("taruwind","bilinear")
  layers <- vector("list",length(timestamps)*length(variables))
  k <- 1L
  for(time_index in seq_along(timestamps)) for(j in seq_along(variables)) {
    time <- timestamps[time_index]
    subset <- wind[wind$timestamp_utc==time,]
    order <- match(sites$site_id,subset$site_id)
    if(anyNA(order)) stop("Area result lacks an output cell")
    layer <- terra::setValues(template,subset[[variables[j]]][order])
    names(layer) <- paste0(prefixes[j],"_",format(time,"%Y%m%dT%H%M%SZ",tz="UTC"))
    layers[[k]] <- layer
    k <- k+1L
  }
  list(wind=wind,maps=do.call(c,layers),
    grid=list(centre=centre,sites=sites,crs=local_crs,
      requested_area_m=area_m,actual_area_m=n*map_resolution,
      map_resolution_m=map_resolution,terrain_resolution_m=resolution,
      terrain_domain_m=domain_m))
}

#' Plot aligned terrain or an hourly wind map
#'
#' Uses base `terra` graphics. Pass the two-layer output of
#' [taruwind_site_layers()] or the `maps` component of
#' [taruwind_downscale_area()]. Layer names are available with `names(x)`.
#' When `layer` is omitted, all supplied layers are drawn.
#' @param x A `terra::SpatRaster` of aligned site layers or wind maps.
#' @param layer One or more layer names or indices; all by default.
#' @param file Optional PNG filename. An existing file is not overwritten.
#' @param width,height PNG dimensions in pixels when `file` is supplied.
#' @return The selected `terra::SpatRaster`, invisibly; draws the plot.
#' @export
taruwind_plot_map <- function(x,layer=NULL,file=NULL,width=1200,height=700) {
  .taruwind_require("terra","map plotting")
  if(!inherits(x,"SpatRaster")) stop("x must be a terra SpatRaster")
  if(!is.null(layer)) x <- x[[layer]]
  if(terra::nlyr(x)<1L) stop("Select at least one map layer")
  if(!is.null(file)) {
    if(!is.character(file)||length(file)!=1L||is.na(file)||!nzchar(file))
      stop("file must be one PNG filename")
    if(file.exists(file)) stop("output already exists: ",file)
    if(!is.numeric(width)||length(width)!=1L||!is.finite(width)||width<=0||
       !is.numeric(height)||length(height)!=1L||!is.finite(height)||height<=0)
      stop("width and height must be positive finite pixel counts")
    grDevices::png(file,width=width,height=height)
    on.exit(grDevices::dev.off(),add=TRUE)
  }
  if(terra::nlyr(x)>1L) {
    old <- graphics::par(mfrow=c(1L,terra::nlyr(x)))
    on.exit(graphics::par(old),add=TRUE,after=FALSE)
  }
  for(j in seq_len(terra::nlyr(x))) {
    name <- names(x)[j]
    pal <- if(name=="elevation_m") grDevices::terrain.colors(80) else
      if(name=="roughness_m") grDevices::hcl.colors(80,"Viridis") else
      grDevices::hcl.colors(80,"YlOrRd",rev=TRUE)
    title <- if(name=="elevation_m") "Elevation (m)" else
      if(name=="roughness_m") "Aerodynamic roughness z0 (m)" else {
        kind <- if(startsWith(name,"taruwind_")) "TaruWind" else "ERA5 bilinear"
        stamp <- sub("^[^_]+_", "", name)
        paste0(kind," speed (m/s), ",substr(stamp,1,4),"-",
          substr(stamp,5,6),"-",substr(stamp,7,8)," ",
          substr(stamp,10,11),":",substr(stamp,12,13)," UTC")
      }
    terra::plot(x[[j]],main=title,col=pal,axes=TRUE,
      xlab="Local east (m)",ylab="Local north (m)")
  }
  invisible(x)
}
