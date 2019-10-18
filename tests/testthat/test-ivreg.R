m.fitivreg <- with(Kmenta, ivreg.fit(cbind(1, P, D), Q, cbind(1, D, F, A)))
m.ivreg <- ivreg(Q ~ P + D | D + F + A, data=Kmenta)

test_that("2SLS coefficients are computed correctly", {
  expect_equal(as.vector(coef(m.fitivreg)), as.vector(coef(m.ivreg)))
})

test_that("2SLS coefficient covariance matrix is computed correctly", {
  expect_equal(as.vector(m.fitivreg$sigma^2*m.fitivreg$cov.unscaled), as.vector(vcov(m.ivreg)))
})
