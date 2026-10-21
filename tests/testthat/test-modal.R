test_that("current modal preparation preserves flat matched-source wind", {
  dem = matrix(0,9,11); rough = matrix(.03,9,11)
  op = taruwind_modal_prepare(dem,rough,90,0,0,450,360,900,.03)
  expect_equal(unname(op$matrix),diag(2),tolerance=1e-12)
  expect_equal(unname(op$surface$factor),rep(1,72),tolerance=1e-12)
  expect_equal(taruwind_modal_predict(op,c(0,3,NA),c(0,4,1))$speed,
               c(0,5,NA),tolerance=1e-12)
  expect_true(all(op$residual==0))
  expect_identical(op$provenance$smoothing_m,0)
})

test_that("current modal terrain response respects datum and scaling invariance", {
  dem = outer(seq(-1,1,length.out=9),seq(-1,1,length.out=11),
               function(y,x) 30*exp(-5*(x*x+y*y)))
  rough = matrix(.03,9,11)
  make = function(h) taruwind_modal_prepare(h,rough,90,0,0,390,330,900,.03)
  op = make(dem); raised = make(dem+2000)
  expect_equal(op$matrix,raised$matrix,tolerance=1e-10)
  u = c(-3,2,0,7);v = c(4,-1,0,-2)
  expect_equal(taruwind_modal_predict(op,2*u,2*v)$speed,
               2*taruwind_modal_predict(op,u,v)$speed,tolerance=1e-12)
  expect_true(all(op$residual<=1.01e-6))
})

test_that("invalid modal input is rejected before native solves", {
  h=matrix(0,5,5);r=matrix(.03,5,5)
  expect_error(taruwind_modal_prepare(h,r,90,0,0,999,180,900,.03),"outside")
  expect_error(taruwind_modal_prepare(h,r,90,0,0,180,180,900,0),"bounds")
  expect_error(taruwind_modal_prepare(h,r,90,0,0,180,180,900,.03,height=101),"bounds")
  expect_error(taruwind_modal_prepare(h,r*NA,90,0,0,180,180,900,.03),"bounds")
})
