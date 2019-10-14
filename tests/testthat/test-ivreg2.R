m.fitivreg2 <- with(Kmenta, fitivreg2(Q, cbind(1, P, D), cbind(1, D, F, A)))
m.ivreg2 <- ivreg2(Q ~ P + D, instruments = ~ D + F + A, data=Kmenta)

test_that("2SLS coefficients are computed correctly", {
  expect_equal(as.vector(coef(m.fitivreg2)), as.vector(coef(m.ivreg2)))
})

test_that("2SLS coefficient covariance matrix is computed correctly", {
  expect_equal(as.vector(m.fitivreg2$vcov), as.vector(vcov(m.ivreg2)))
})
