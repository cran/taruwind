test_that("site input maps use the same aligned grid as the solver", {
  skip_if_not_installed("terra")
  crs = "+proj=aeqd +lat_0=41 +lon_0=29 +datum=WGS84 +units=m"
  dem = terra::rast(nrows=41,ncols=41,xmin=-1845,xmax=1845,
                    ymin=-1845,ymax=1845,crs=crs)
  terra::values(dem) = 0
  rough = dem*0 + .03
  layers = taruwind_site_layers(29,41,roughness=rough,dem=dem,domain_m=450)
  expect_identical(names(layers),c("elevation_m","roughness_m"))
  expect_equal(terra::nrow(layers),5)
  expect_equal(terra::values(layers[["roughness_m"]],mat=FALSE),rep(.03,25),tolerance=1e-7)
  file = tempfile(fileext=".png")
  taruwind_plot_map(layers,layer="elevation_m",file=file)
  expect_true(file.exists(file))
  expect_error(taruwind_plot_map(layers,layer="elevation_m",file=file),"exists")
  both_file = tempfile(fileext=".png")
  taruwind_plot_map(layers,file=both_file)
  expect_true(file.exists(both_file))
  op = taruwind_prepare_site(29,41,rough,900,.03,dem=dem,domain_m=450)
  expect_equal(taruwind_predict(op,3,4)$speed,5,tolerance=1e-6)
})

test_that("area output keeps each target and hour in its raster cell", {
  skip_if_not_installed("terra")
  crs = "+proj=aeqd +lat_0=41 +lon_0=29 +datum=WGS84 +units=m"
  dem = terra::rast(nrows=41,ncols=41,xmin=-1845,xmax=1845,
                    ymin=-1845,ymax=1845,crs=crs)
  terra::values(dem) = 0
  rough = dem*0 + .03
  local_mocked_bindings(.taruwind_era5_axis=function(...) NULL,
    .taruwind_download_era5_ready=function(target,period,...) {
    time = seq(period$start,period$end,by="hour")
    data = target[rep(seq_len(nrow(target)),each=length(time)),]
    data$timestamp_utc = rep(time,nrow(target))
    data$u = seq_len(nrow(data))
    data$v = 0
    data$wind_speed = data$u
    tibble::as_tibble(data)
  })
  result = taruwind_downscale_area(29,41,
    start="2020-01-01 00:00:00",end="2020-01-01 01:00:00",
    area_m=990,map_resolution=500,dem=dem,roughness=rough,
    source_z0=.03,depth=900,domain_m=450,workers=1,quiet=TRUE)
  expect_equal(nrow(result$wind),18)
  expect_equal(terra::nlyr(result$maps),4)
  expect_equal(result$grid$actual_area_m,1500)
  expect_equal(terra::values(result$maps[["bilinear_20200101T000000Z"]],mat=FALSE),
               seq(1,17,by=2))
  expect_equal(terra::values(result$maps[["taruwind_20200101T000000Z"]],mat=FALSE),
               seq(1,17,by=2),tolerance=1e-6)
  reused = taruwind_downscale_area(29,41,
    start="2020-01-02 00:00:00",end="2020-01-02 01:00:00",
    area_m=990,map_resolution=500,responses=attr(result$wind,"responses"),
    domain_m=450,quiet=TRUE)
  expect_equal(reused$wind$taruwind_speed,result$wind$taruwind_speed)
  expect_true("taruwind_20200102T000000Z" %in% names(reused$maps))
  expect_error(taruwind_downscale_area(29,41,area_m=2000,map_resolution=100),
               "Output map has")
  expect_error(taruwind_downscale_area(29,41,
    period=c("2020-01-01","2020-01-08"),area_m=990,map_resolution=500),
    "168 hourly timestamps")
})
