test_that("UTC dates retain inclusive hourly boundaries", {
  p=.taruwind_resolve_period(period=c("2020-01-01","2020-01-01"))
  expect_equal(as.numeric(difftime(p$end,p$start,units="hours")),23)
  expect_error(.taruwind_resolve_period(start="2020-01-01",period=c("2020-01-01","2020-01-01")),"combined")
  expect_error(.taruwind_period("2020-01-01 00:15:00","2020-01-01 01:00:00"),"whole hours")
  expect_error(.taruwind_period("1939-12-31","1940-01-01"),"before 1940")
  expect_error(.taruwind_period("2999-01-01","2999-01-02"),"future")
  expect_equal(.taruwind_request_hours(p),24)
  expect_equal(.taruwind_request_hours(
    .taruwind_period("2020-01-01","2020-01-07")),
    taruwind_limits()$max_area_hours)
})

test_that("dynamic ERA5 availability checks every requested hour", {
  hours = seq(as.POSIXct("2020-01-01 00:00:00",tz="UTC"),
    by="hour",length.out=4)
  p = .taruwind_period("2020-01-01 00:00:00","2020-01-01 02:00:00")
  expect_equal(.taruwind_check_era5_period(p,hours),1:3)
  expect_error(.taruwind_check_era5_period(
    .taruwind_period("2020-01-01 02:00:00","2020-01-01 05:00:00"),
    hours),"2020-01-01 03:00 UTC")
  expect_error(.taruwind_check_era5_period(p,hours[-2]),"unavailable")
  expect_error(.taruwind_check_era5_period(p,as.POSIXct(character(),tz="UTC")),
    "empty or invalid")
})

test_that("request guard counts reflected terrain and rejects excessive networks", {
  p = .taruwind_period("2020-01-01 00:00:00","2020-01-01 00:00:00")
  small = data.frame(site_id=paste0("s",seq_len(59)))
  large = data.frame(site_id=paste0("s",seq_len(60)))
  expect_silent(.taruwind_guard_request(small,p,resolution=90,domain_km=30))
  expect_error(.taruwind_guard_request(large,p,resolution=90,domain_km=30,
    allow_override=FALSE),"total reflected computational cells")
  expect_error(.taruwind_guard_request(small,p,resolution=0,domain_km=30),
    "positive finite")
})

test_that("automatic requests are bounded by ten calendar years", {
  site = data.frame(site_id="one")
  ten_years = .taruwind_period("2016-01-01","2025-12-31")
  too_long = .taruwind_period("2016-01-01","2026-01-01")
  expect_silent(.taruwind_guard_request(site,ten_years))
  expect_error(.taruwind_guard_request(site,too_long,allow_override=FALSE),
    "10 calendar years")
  expect_warning(.taruwind_guard_request(site,too_long,allow_large=TRUE),
    "10 calendar years")
})

test_that("the main workflow checks ERA5 availability before terrain or profiles", {
  local_mocked_bindings(
    .taruwind_era5_axis=function(...) stop("not available in ERA5"),
    taruwind_source_profile=function(...) stop("unexpected source-profile request")
  )
  expect_error(taruwind_downscale(lon=29,lat=41,
    period=c("2020-01-01","2020-01-01"),roughness="worldcover",
    quiet=TRUE),"not available in ERA5")
})

test_that("local timestamps do not change UTC instants", {
  skip_if_not_installed("lutz")
  sites=.taruwind_sites(lon=29,lat=41)
  x=data.frame(site_id=sites$site_id,timestamp_utc=as.POSIXct("2020-01-01",tz="UTC"))
  y=.taruwind_add_local_time(x,sites,TRUE)
  expect_identical(y$timestamp_utc,x$timestamp_utc)
  expect_equal(y$time_zone,"Europe/Istanbul")
  expect_equal(y$timestamp_local,"2020-01-01T03:00:00+0300")
})
