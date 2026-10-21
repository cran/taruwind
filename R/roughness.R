#' Download WorldCover roughness on a target raster grid
#'
#' Reads windows from public ESA WorldCover 2021 v200 COGs, samples categorical
#' land cover by nearest neighbour, and applies the WindKit 2.1.0 default
#' WorldCover lookup. This is an empirical land-cover conversion, not a direct
#' roughness measurement. It is not interchangeable with the manuscript's
#' CORINE-based inputs. Water is floored at 0.0002 m for logarithmic profiles.
#' @param template Single-layer projected SpatRaster defining the target grid.
#' @param cache_dir Directory for cached categorical windows.
#' @return Roughness SpatRaster in metres, with conversion provenance.
#' @references ESA WorldCover: <https://esa-worldcover.org/en/data-access>.
#'   WindKit roughness tables:
#'   <https://docs.wasp.dk/windkit/latest/user_guide/topography/roughness_maps.html>.
#' @export
taruwind_download_worldcover <- function(template,
    cache_dir=file.path(file.path(tempdir(), "taruwind-cache"),"worldcover")) {
  .taruwind_require("terra","WorldCover window extraction")
  if(!inherits(template,"SpatRaster")||terra::nlyr(template)!=1L||
     !nzchar(terra::crs(template))||terra::is.lonlat(template)||terra::ncell(template)>2e6)
    stop("Projected single-layer template with at most two million cells required")
  dir.create(cache_dir,recursive=TRUE,showWarnings=FALSE)
  ll <- terra::crds(terra::project(terra::vect(terra::xyFromCell(template,
    seq_len(terra::ncell(template))),crs=terra::crs(template)),"EPSG:4326"))
  owners <- cbind(3*floor(ll[,1]/3),3*floor(ll[,2]/3))
  tiles <- unique(owners)
  if(nrow(tiles)>9) stop("WorldCover request spans too many source tiles")
  classes <- rep(NA_real_,nrow(ll));urls<-character()
  for(i in seq_len(nrow(tiles))) {
    lon<-tiles[i,1];lat<-tiles[i,2]
    code<-sprintf("%s%02d%s%03d",if(lat>=0) "N" else "S",abs(lat),if(lon>=0) "E" else "W",abs(lon))
    url<-paste0("https://esa-worldcover.s3.eu-central-1.amazonaws.com/v200/2021/map/ESA_WorldCover_10m_2021_v200_",code,"_Map.tif")
    # Grid signature prevents reuse on different projections or alignments.
    sig<-list(crs=terra::crs(template),extent=as.vector(terra::ext(template)),dim=dim(template))
    # Each grid receives its own MD5-keyed subdirectory.
    tmp<-tempfile(tmpdir=cache_dir);saveRDS(sig,tmp)
    key<-unname(tools::md5sum(tmp));unlink(tmp)
    folder<-file.path(cache_dir,key);dir.create(folder,showWarnings=FALSE)
    path<-file.path(folder,paste0(code,".tif"))
    if(!file.exists(path)) {
      src<-terra::rast(paste0("/vsicurl/",url))
      sampled<-terra::project(src,template,method="near",threads=FALSE)
      terra::writeRaster(sampled,path,overwrite=FALSE,gdal="COMPRESS=DEFLATE")
      saveRDS(sig,file.path(folder,"grid_signature.rds"))
    }
    cached<-terra::rast(path)
    if(!terra::compareGeom(cached,template,stopOnError=FALSE)) stop("Cached WorldCover geometry mismatch")
    owned<-owners[,1]==lon & owners[,2]==lat
    classes[owned]<-terra::values(cached,mat=FALSE)[owned]
    urls<-c(urls,url)
  }
  lookup<-c(`10`=1.5,`20`=.3,`30`=.03,`40`=.1,`50`=.5,`60`=.005,
            `70`=.003,`80`=.0002,`90`=.1,`95`=.7,`100`=.01)
  if(anyNA(classes)||any(!classes%in%as.numeric(names(lookup))))
    stop("Incomplete or unknown WorldCover classes; no silent fallback is applied")
  result<-terra::setValues(template,unname(lookup[as.character(classes)]))
  attr(result,"taruwind_roughness_provenance")<-list(product="WorldCover 2021 v200",
    conversion="WindKit 2.1.0 default WorldCover lookup",water_floor_m=.0002,
    resampling="nearest category",urls=urls)
  result
}
