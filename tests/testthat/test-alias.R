## data with collinear z1 and z2
set.seed(0)
d <- data.frame(
  y = rnorm(20),
  x1 = sin(1:20),
  x2 = cos(1:20),
  z1 = 1:20,
  z2 = 2:21/4
)

test_that("collinear regressors", {
  m0 <- expect_silent(ivreg(y ~ z1 | x1 + x2, data = d))
  m <- expect_silent(ivreg(y ~ z1 + z2 | x1 + x2, data = d))
  expect_true(summary(m)$aliased["z2"])
  expect_equal(coef(m0), coef(m)[!is.na(coef(m))])
  expect_equal(residuals(m0), residuals(m))
})

test_that("more regressors than instruments", {
  m0 <- expect_silent(ivreg(y ~ x1 | z1, data = d))
  m <- expect_warning(ivreg(y ~ x1 + x2 | z1, data = d))
  expect_true(summary(m)$aliased["x2"])
  expect_equal(coef(m0), coef(m)[!is.na(coef(m))])
  expect_equal(residuals(m0), residuals(m))
})

test_that("collinear instrumental variables", {
  m0 <- expect_silent(ivreg(y ~ x1 | z1, data = d))
  m <- expect_warning(ivreg(y ~ x1 | z1 + z2, data = d))
  expect_true(all(!summary(m, diagnostics = FALSE)$aliased))
  expect_equal(coef(m0), coef(m))
  expect_equal(residuals(m0), residuals(m))
})

test_that("collinear instrumental variables 2", {
  m0 <- expect_silent(ivreg(y ~ x1 | z1, data = d))
  m <- expect_warning(ivreg(y ~ x1 + x2 | z1 + z2, data = d))
  expect_true(summary(m, diagnostics = FALSE)$aliased["x2"])
  expect_equal(coef(m0), coef(m)[!is.na(coef(m))])
  expect_equal(residuals(m0), residuals(m))
})
