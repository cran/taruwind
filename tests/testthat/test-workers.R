test_that("PSOCK workers use the requested package installation", {
  h=matrix(0,5,5)
  q=list(dem=h,roughness=h+.03,dx=90,x0=0,y0=0,x=180,y=180,depth=900,source_z0=.03)
  answers=taruwind_build_many(list(q,q),workers=2)
  expect_length(answers,2)
  expect_identical(answers[[1]]$provenance$package_version,"0.1.1")
  expect_equal(answers[[1]]$matrix,answers[[2]]$matrix,tolerance=1e-12)
})
