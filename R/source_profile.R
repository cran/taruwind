#' Compute the model's neutral modal depth from latitude
#' @param lat Numeric latitude(s) in degrees. Automatic use is restricted to
#'   absolute latitudes 5 through 85 degrees; provide explicit depth elsewhere.
#' @return Neutral depth in metres, using the study convention
#'   `0.2 * (0.4 * 10 / log(10 / 0.01)) / abs(1.4544e-4 * sin(latitude))`.
#'   This is a fixed numerical scale, not measured or hourly ERA5 BLH.
#' @export
taruwind_neutral_depth <- function(lat) {
  if(!is.numeric(lat)||!length(lat)||any(!is.finite(lat))||any(abs(lat)<5|abs(lat)>85))
    stop("Automatic neutral depth requires absolute latitude 5 to 85 degrees; supply depth explicitly elsewhere")
  .2*(.4*10/log(10/.01))/abs(1.4544e-4*sin(lat*pi/180))
}

#' Obtain the explicit ERA5 source-profile convention
#'
#' Retrieves the four-node ERA5 forecast surface roughness (FSR) climatology
#' from CDS monthly means, averages the 12 months and 24 hours equally and
#' bilinearly interpolates the result. The default reference year is 2021,
#' matching the study input convention, not necessarily the requested wind year.
#' At 10 m the effective diagnostic reference is bounded to [0.0002, 0.03] m.
#' This cap is a source-profile convention, not a claim that ERA5's actual
#' roughness is globally 0.03 m. At 100 m the raw FSR is retained with a 0.0002 m
#' floor. The surface multiplier is unity at the 100 m matching reference.
#'
#' CDS credentials and dataset terms acceptance are required. A queued job is
#' cached and resumed on the next call; timeout does not submit duplicate jobs.
#' No Python, keyring or observations are needed. Existing FSR NetCDF files may
#' be supplied for offline use. REST access is subject to CDS service changes.
#' @param lon,lat One target coordinate in degrees.
#' @param wind_height 10 or 100 metres.
#' @param reference_year FSR climatology year, default 2021.
#' @param fsr_file Optional existing NetCDF containing `fsr`, longitude,
#'   latitude and optional time dimensions. Its actual averaging period is the
#'   caller's responsibility and is not inferred from `reference_year`.
#' @param api_key CDS key; otherwise standard environment/config discovery.
#' @param cache_dir Local request, job and NetCDF cache directory.
#' @param wait_seconds Maximum queue polling time per call, default 60 seconds.
#' @param quiet Suppress cache and CDS queue messages.
#' @return A list containing raw/effective roughness, neutral depth and provenance.
#' @references <https://cds.climate.copernicus.eu/how-to-api>
#' @export
taruwind_source_profile <- function(lon,lat,wind_height=10,reference_year=2021,
    fsr_file=NULL,api_key=NULL,
    cache_dir=file.path(file.path(tempdir(), "taruwind-cache"),"source_profile"),
    wait_seconds=60,quiet=TRUE) {
  .taruwind_require("ncdf4","ERA5 source-profile reading")
  site<-.taruwind_sites(lon=lon,lat=lat)
  if(nrow(site)!=1L || lon>=179.75 || lon< -179.75) stop("One non-antimeridian target required")
  height<-.taruwind_wind_height(wind_height)
  depth<-taruwind_neutral_depth(lat)
  if(!is.numeric(reference_year)||length(reference_year)!=1||!is.finite(reference_year)||
     reference_year<1940||reference_year!=floor(reference_year)) stop("Valid integer reference_year required")
  supplied<-!is.null(fsr_file)
  if(!supplied) {
    if(length(wait_seconds)!=1||!is.finite(wait_seconds)||wait_seconds<0) stop("Non-negative wait_seconds required")
    west<-floor(lon*4)/4;south<-floor(lat*4)/4
    request<-list(product_type=list("monthly_averaged_reanalysis_by_hour_of_day"),
      variable=list("forecast_surface_roughness"),year=list(as.character(reference_year)),
      month=as.list(sprintf("%02d",1:12)),time=as.list(sprintf("%02d:00",0:23)),
      area=c(south+.25,west,south,west+.25),data_format="netcdf",download_format="unarchived")
    folder<-file.path(cache_dir,sprintf("fsr_%d_%+.2f_%+.2f",reference_year,west,south))
    dir.create(folder,recursive=TRUE,showWarnings=FALSE)
    fsr_file<-file.path(folder,"fsr.nc")
    if (!quiet) message("  FSR ", reference_year, ": ", if (file.exists(fsr_file))
      "using cached NetCDF" else "requesting monthly-by-hour reference data from CDS")
    if(!file.exists(fsr_file)) .taruwind_cds_fsr(request,folder,api_key,wait_seconds,quiet)
  }
  nc<-ncdf4::nc_open(fsr_file);on.exit(ncdf4::nc_close(nc),add=TRUE)
  if(!"fsr"%in%names(nc$var)) stop("NetCDF must contain the fsr variable")
  dims<-vapply(nc$var$fsr$dim,`[[`,"","name")
  ilon<-match("longitude",dims);ilat<-match("latitude",dims)
  if(anyNA(c(ilon,ilat))) stop("FSR requires longitude and latitude dimensions")
  a<-ncdf4::ncvar_get(nc,"fsr",collapse_degen=FALSE)
  if(any(!is.finite(a))||any(a<=0)) stop("FSR contains invalid values")
  if(!nc$var$fsr$units%in%c("m","metres","meters")) stop("FSR units must be metres")
  m<-apply(a,c(ilon,ilat),mean)
  lons<-ncdf4::ncvar_get(nc,"longitude")%%360;lats<-ncdf4::ncvar_get(nc,"latitude")
  axis_weights<-function(values,target) {
    exact<-which(abs(values-target)<1e-8)
    if(length(exact)) return(list(index=exact[1],weight=1))
    below<-which(values<target);above<-which(values>target)
    if(!length(below)||!length(above)) stop("Target outside FSR coordinate coverage")
    lo<-below[which.max(values[below])];hi<-above[which.min(values[above])]
    w<-(target-values[lo])/(values[hi]-values[lo])
    list(index=c(lo,hi),weight=c(1-w,w))
  }
  ax<-axis_weights((lons-lon+180)%%360-180+lon,lon);ay<-axis_weights(lats,lat)
  raw<-sum(m[ax$index,ay$index,drop=FALSE]*outer(ax$weight,ay$weight))
  effective<-if(height==10) max(.0002,min(.03,raw)) else max(.0002,raw)
  list(depth=depth,source_z0=effective,raw_fsr=raw,
    provenance=list(adapter="ERA5 diagnostic reference",wind_height=height,
      reference_year=if(supplied) NA_integer_ else reference_year,
      averaging="equal-weight mean of supplied FSR time slices; bilinear four-node interpolation",
      convention=if(height==10) "FSR bounded to [0.0002, 0.03] m" else "FSR floored at 0.0002 m",
      file=normalizePath(fsr_file),md5=unname(tools::md5sum(fsr_file)),observations_used=FALSE))
}

.taruwind_cds_fsr <- function(request,folder,api_key,wait_seconds,quiet=TRUE) {
  key<-.taruwind_cds_key(api_key)
  base<-"https://cds.climate.copernicus.eu/api/retrieve/v1"
  call<-function(path,body=NULL) {
    h<-curl::new_handle(connecttimeout=20,timeout=60)
    curl::handle_setheaders(h,`PRIVATE-TOKEN`=key,`Content-Type`="application/json")
    if(!is.null(body)) curl::handle_setopt(h,postfields=jsonlite::toJSON(body,auto_unbox=TRUE))
    r<-curl::curl_fetch_memory(paste0(base,path),h)
    if(r$status_code>=400) stop("CDS source-profile request failed (HTTP ",r$status_code,
      "); check dataset terms, credentials and service status")
    jsonlite::fromJSON(rawToChar(r$content),simplifyVector=FALSE)
  }
  jf<-file.path(folder,"job.json")
  if(file.exists(jf)) job<-jsonlite::fromJSON(jf) else {
    job<-call("/processes/reanalysis-era5-single-levels-monthly-means/execution",list(inputs=request))
    if(is.null(job$jobID)) stop("CDS response did not include jobID")
    jsonlite::write_json(list(jobID=job$jobID),jf,auto_unbox=TRUE)
    jsonlite::write_json(request,file.path(folder,"request.json"),auto_unbox=TRUE,pretty=TRUE)
  }
  begin<-Sys.time()
  previous_status = NULL
  repeat {
    s<-call(paste0("/jobs/",job$jobID))
    if (!quiet && !identical(previous_status,s$status)) message("  CDS source-profile job: ",s$status)
    previous_status = s$status
    if(identical(s$status,"successful")) break
    if(s$status%in%c("failed","dismissed")) stop("Cached CDS source-profile job ",s$status,"; inspect ",jf)
    if(as.numeric(difftime(Sys.time(),begin,units="secs"))>=wait_seconds)
      stop("CDS source-profile job is still queued/running. Repeat this call to resume the cached job: ",jf)
    Sys.sleep(min(5,max(.1,wait_seconds)))
  }
  result<-call(paste0("/jobs/",job$jobID,"/results"))
  url<-result$asset$value$href
  if(is.null(url)||!grepl("^https://",url)) stop("CDS returned no HTTPS result asset")
  tmp<-file.path(folder,"fsr.partial.nc")
  curl::curl_download(url,tmp,quiet=TRUE,handle=curl::new_handle(connecttimeout=20,timeout=120))
  nc<-ncdf4::nc_open(tmp);ncdf4::nc_close(nc)
  if(!file.rename(tmp,file.path(folder,"fsr.nc"))) stop("Could not finalize FSR cache")
  invisible(NULL)
}
