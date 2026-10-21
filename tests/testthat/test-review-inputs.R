test_that("factor coordinates use labels, never internal integer codes", {
  s = .taruwind_sites(data.frame(lon=factor(c("29.02","28.9")),
                                lat=factor(c("41.02","40.1"))))
  expect_equal(s$lon,c(29.02,28.9))
  expect_equal(s$lat,c(41.02,40.1))
  expect_error(.taruwind_sites(data.frame(lon=factor("east"),lat=41)),"finite coordinates")
})

test_that("provided site IDs are nonempty and nonmissing", {
  for (id in list(NA_character_,"","  ")) {
    expect_error(.taruwind_sites(data.frame(lon=29,lat=41,site_id=id)),"site identifiers")
  }
  expect_equal(.taruwind_sites(lon=29,lat=41)$site_id,"site_1")
})

test_that("invalid local timezone cannot silently label UTC as local time", {
  sites = .taruwind_sites(lon=29,lat=41)
  x = data.frame(site_id=sites$site_id,timestamp_utc=as.POSIXct("2020-01-01",tz="UTC"))
  expect_error(.taruwind_add_local_time(x,sites,"Not/A_Timezone"),"valid Olson")
  expect_error(.taruwind_add_local_time(x,sites,NA_character_),"local_tz")
  y = .taruwind_add_local_time(x,sites,"Europe/Istanbul")
  expect_equal(y$timestamp_local,"2020-01-01T03:00:00+0300")
  expect_identical(y$timestamp_utc,x$timestamp_utc)
})
