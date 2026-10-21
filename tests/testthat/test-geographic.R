test_that("geographic preparation preserves flat uniform flow", {
  skip_if_not_installed("terra")
  r = terra::rast(nrows=21,ncols=21,xmin=28.9,xmax=29.1,ymin=40.9,ymax=41.1,crs="EPSG:4326")
  terra::values(r) = 0
  op = taruwind_prepare_site(29,41,r+.03,900,.03,dem=r,domain_m=450)
  expect_equal(taruwind_predict(op,3,4)$speed,5,tolerance=1e-6)
  expect_equal(op$provenance$geographic$actual_domain_m,450)
  expect_error(taruwind_prepare_site(29,41,r+.03,900,.03,dem=r,domain_m=1e7),"limit")
  expect_error(taruwind_downscale(lon=29,lat=41,period=c("2020-01-01","2020-01-01")),"roughness")
})

test_that("raster row orientation agrees with explicit metric-grid preparation", {
  skip_if_not_installed("terra")
  crs="+proj=aeqd +lat_0=41 +lon_0=29 +datum=WGS84 +units=m"
  r=terra::rast(nrows=7,ncols=7,xmin=-315,xmax=315,ymin=-315,ymax=315,crs=crs)
  xy=terra::xyFromCell(r,1:49)
  terra::values(r)=20*exp(-((xy[,1]-90)^2+(xy[,2]-90)^2)/20000)
  z=r*0+.03
  a=taruwind_prepare_site(29,41,z,900,.03,dem=r,domain_m=630)
  m=terra::as.matrix(r,wide=TRUE)
  b=taruwind_build(m[7:1,,drop=FALSE],matrix(.03,7,7),90,-270,-270,0,0,900,.03)
  expect_equal(a$matrix,b$matrix,tolerance=1e-9)
  expect_equal(taruwind_predict(a,c(3,-3),c(4,-4)),
    taruwind_predict(b,c(3,-3),c(4,-4)),tolerance=1e-9)
})

test_that("download guards reject broad or invalid input without HTTP", {
  skip_if_not_installed("terra")
  expect_error(taruwind_download_copernicus_dem(c(-181,0,1,1),"EPSG:32635"),"bbox")
  expect_error(taruwind_download_copernicus_dem(c(0,0,40,40),"EPSG:32635"),"16 DEM")
  expect_error(taruwind_download_worldcover(terra::rast()),"Projected")
})

test_that("geographic pipeline joins forcing to the current response", {
  skip_if_not_installed("terra")
  r = terra::rast(nrows=21,ncols=21,xmin=28.9,xmax=29.1,ymin=40.9,ymax=41.1,crs="EPSG:4326")
  terra::values(r)=0
  local_mocked_bindings(.taruwind_era5_axis=function(...) NULL,
    .taruwind_download_era5_ready=function(target,...) {
    tibble::tibble(site_id=target$site_id,lon=target$lon,lat=target$lat,
      timestamp_utc=as.POSIXct("2020-01-01",tz="UTC"),u=3,v=4,wind_speed=5)
  },taruwind_source_profile=function(...) list(source_z0=.03,
    provenance=list(adapter="test fixture")))
  x = taruwind_downscale(lon=29,lat=41,period=c("2020-01-01","2020-01-01"),
    dem=r,roughness=r+.03,depth=900,source_z0=.03,domain_m=450)
  expect_equal(x$taruwind_speed,5,tolerance=1e-6)
  y = taruwind_downscale(sites=matrix(c(29,41),1,dimnames=list(NULL,c("lon","lat"))),
    period=c("2020-01-01","2020-01-01"),responses=attr(x,"responses"))
  expect_equal(y$taruwind_speed,x$taruwind_speed)
  expect_error(taruwind_downscale(lon=30,lat=41,period=c("2020-01-01","2020-01-01"),
    responses=attr(x,"responses")),"coordinates")
  pair = taruwind_downscale(lon=c(29,29.001),lat=c(41,41.001),
    period=c("2020-01-01","2020-01-01"),dem=r,roughness=r+.03,
    depth=900,source_z0=.03,domain_m=450,workers=2)
  expect_equal(pair$taruwind_speed,c(5,5),tolerance=1e-6)
  expect_equal(names(attr(pair,"responses")),c("site_1","site_2"))
  auto = taruwind_downscale(lon=29,lat=41,
    period=c("2020-01-01","2020-01-01"),dem=r,roughness=r+.03,domain_m=450)
  expect_equal(auto$taruwind_speed,5,tolerance=1e-6)
  expect_equal(attr(auto,"responses")[[1]]$provenance$depth,taruwind_neutral_depth(41))
  expect_equal(attr(auto,"responses")[[1]]$provenance$source_profile$adapter,"test fixture")
})
