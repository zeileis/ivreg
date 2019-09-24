m.fitivreg <- with(Kmenta, fitivreg(Q, cbind(1, P, D), cbind(1, D, F, A)))
m.ivreg <- ivreg(Q ~ P + D, instruments = ~ D + F + A, data=Kmenta)

test_that("2SLS coefficients are computed correctly", {
  expect_equal(as.vector(coef(m.fitivreg)), as.vector(coef(m.ivreg)))
})

test_that("2SLS coefficient covariance matrix is computed correctly", {
  expect_equal(as.vector(m.fitivreg$vcov), as.vector(vcov(m.ivreg)))
})
