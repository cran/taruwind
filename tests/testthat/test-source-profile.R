test_that("neutral depth matches the study formula and refuses equatorial singularity", {
  expect_equal(taruwind_neutral_depth(41),.2*(4/log(1000))/(1.4544e-4*sin(41*pi/180)))
  expect_equal(taruwind_neutral_depth(-41),taruwind_neutral_depth(41))
  expect_error(taruwind_neutral_depth(0),"latitude")
})
test_that("FSR averaging, spatial interpolation and diagnostic convention are explicit", {
  skip_if_not_installed("ncdf4")
  f=tempfile(fileext=".nc");on.exit(unlink(f),add=TRUE)
  dims=list(ncdf4::ncdim_def("longitude","degrees_east",c(29,29.25)),
    ncdf4::ncdim_def("latitude","degrees_north",c(41.25,41)),
    ncdf4::ncdim_def("time","hours since 2021-01-01",0:1))
  nc=ncdf4::nc_create(f,ncdf4::ncvar_def("fsr","m",dims,prec="double"))
  a=array(c(.01,.05,.09,.13,.03,.07,.11,.15),c(2,2,2))
  ncdf4::ncvar_put(nc,"fsr",a);ncdf4::nc_close(nc)
  p=taruwind_source_profile(29.125,41.125,fsr_file=f)
  expect_equal(p$raw_fsr,mean(a))
  expect_equal(p$source_z0,.03)
  expect_true(is.na(p$provenance$reference_year))
  p100=taruwind_source_profile(29.125,41.125,wind_height=100,fsr_file=f)
  expect_equal(p100$source_z0,mean(a))
  expect_equal(taruwind_source_profile(29,41.25,fsr_file=f)$raw_fsr,.02)
})
