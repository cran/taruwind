test_that("cache inspection is read-only and includes hidden metadata", {
  root = tempfile(); dir.create(root)
  on.exit(unlink(root,recursive=TRUE),add=TRUE)
  writeBin(as.raw(1:3),file.path(root,".zmetadata"))
  info = taruwind_cache_info(root)
  expect_equal(info$file,".zmetadata")
  expect_equal(attr(info,"total_bytes"),3)
  expect_true(file.exists(file.path(root,".zmetadata")))
  expect_equal(nrow(taruwind_cache_info(file.path(root,"absent"))),0)
})

test_that("quiet workflow preserves output and response reuse", {
  axis_calls = 0L
  local_mocked_bindings(.taruwind_era5_axis=function(...) {
    axis_calls <<- axis_calls + 1L
    list(metadata_refresh=TRUE)
  },
    .taruwind_download_era5_ready=function(
        target,period,axis,local_tz,include_nodes,quiet,wind_height) {
    expect_true(quiet)
    expect_true(axis$metadata_refresh)
    tibble::tibble(site_id=target$site_id,u=3,v=4,wind_speed=5)
  })
  h = matrix(0,5,5)
  op = taruwind_build(h,h+.03,90,0,0,180,180,900,.03)
  op$provenance$geographic = list(lon=29,lat=41)
  responses = list(site_1=op)
  expect_message({x = taruwind_downscale(lon=29,lat=41,
    period=c("2020-01-01","2020-01-01"),responses=responses,quiet=TRUE)},NA)
  expect_equal(axis_calls,1L)
  expect_equal(x$taruwind_speed,5)
})
