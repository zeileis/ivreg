
m.ivreg2 <- with(Kmenta, fitivreg2(Q, cbind(1, P, D), cbind(1, D, F, A)))
m.stage1 <- lm(cbind(1, P, D) ~ cbind(1, D, F, A) - 1, data=Kmenta)
m.stage2 <- lm(Kmenta$Q ~ fitted(m.stage1) - 1)


test_that("2SLS coefficients are computed correctly", {
  expect_equal(as.vector(coef(m.ivreg2)), as.vector(coef(m.stage2)))
})

se.ivreg2 <- as.vector(sqrt(diag(m.ivreg2$vcov)))
se <- as.vector(sqrt(diag(vcov(m.stage2))))
residuals <- as.vector(with(Kmenta, Q - cbind(1, P, D) %*% coef(m.stage2)))
se <- sqrt(sum(residuals^2)/(nrow(Kmenta) - 3))*se/sigma(m.stage2)

test_that("2SLS coefficient standard errors are computed correctly", {
  expect_equal(se.ivreg2, se)
})

test_that("model residuals are computed correctly", {
  expect_equal(residuals, residuals(m.ivreg2))
})


m.ivreg2.w <- with(Kmenta, fitivreg2(Q, cbind(1, P, D), cbind(1, D, F, A), wt=Q))
m.stage1.w <- lm(cbind(1, P, D) ~ cbind(1, D, F, A) - 1, weights=Q, data=Kmenta)
m.stage2.w <- lm(Kmenta$Q ~ fitted(m.stage1.w) - 1, weights=Kmenta$Q)

test_that("2SLS coefficients are computed correctly with weights", {
  expect_equal(as.vector(coef(m.ivreg2.w)), as.vector(coef(m.stage2.w)))
})

se.ivreg2.w <- as.vector(sqrt(diag(m.ivreg2.w$vcov)))
se.w <- as.vector(sqrt(diag(vcov(m.stage2.w))))
residuals <- as.vector(with(Kmenta, Q - cbind(1, P, D) %*% coef(m.stage2.w)))
se.w <- sqrt(sum(Kmenta$Q*residuals^2)/(nrow(Kmenta) - 3))*se.w/sigma(m.stage2.w)

test_that("2SLS coefficient standard errors are computed correctly with weights", {
  expect_equal(se.ivreg2.w, se.w)
})
