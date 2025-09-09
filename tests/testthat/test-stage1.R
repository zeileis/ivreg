set.seed(123)
z1 <- rnorm(100)
z2 <- rnorm(100)
z3 <- rnorm(100)
x1 <- z1 + z2 + rnorm(100)
x2 <- z2 + z3 + rnorm(100)
y <- 10 + z1 + x1 + x2 + rnorm(100)

D <- data.frame(z1, z2, z3, x1, x2, y, const=1)

m1 <- ivreg(y ~ z1 + x1 + x2 | z1 + z2 + z3, data=D)
ml1 <- lm(cbind(x1, x2) ~ z1 + z2 + z3, data=D)

test_that("stage1 coefficients are computed correctly", {
  expect_equal(as.vector(coef(m1, component="stage1")), as.vector(coef(ml1)))
})

test_that("stage1 coefficient covariances are computed correctly", {
  expect_equal(vcov(m1, component="stage1"), vcov(ml1))
})

m2 <- expect_warning(ivreg(y ~ z1 + x1 + x2 + const | z1 + z2 + z3 + const, data=D))
ml2 <- lm(cbind(x1, x2) ~ z1 + z2 + z3 + const, data=D)

test_that("stage1 rank-deficient coefficients are computed correctly", {
  expect_equal(as.vector(coef(m2, component="stage1")), as.vector(coef(ml2)))
})

test_that("stage1 rank-deficient coefficient covariances are computed correctly", {
  expect_equal(vcov(m2, component="stage1"), vcov(ml2))
})

test_that("stage1 rank-deficient coefficients are computed correctly, complete=FALSE", {
  expect_equal(as.vector(coef(m2, component="stage1", complete=FALSE)), 
               as.vector(coef(ml2, complete=FALSE)))
})

test_that("stage1 rank-deficient coefficient covariances are computed correctly, rank=FALSE", {
  expect_equal(vcov(m2, component="stage1", complete=FALSE), 
               vcov(ml2, complete=FALSE))
})


y2 <- 10 + z1 + x1 + rnorm(100)

D2 <- data.frame(z1, z2, z3, x1, y, const=1)

m3 <- ivreg(y2 ~ z1 + x1 | z1 + z2 + z3, data=D2)
ml3 <- lm(x1 ~ z1 + z2 + z3, data=D)

test_that("stage1 coefficients are computed correctly with one endog x", {
  expect_equal(as.vector(coef(m3, component="stage1")), as.vector(coef(ml3)))
})

test_that("stage1 coefficient covariances are computed correctly with one endog x", {
  expect_equal(vcov(m3, component="stage1"), vcov(ml3))
})

m4 <- expect_warning(ivreg(y2 ~ z1 + x1 + const | z1 + z2 + z3 + const, data=D))
ml4 <- lm(x1 ~ z1 + z2 + z3 + const, data=D)

test_that("stage1 rank-deficient coefficients are computed correctly with one endog x", {
  expect_equal(as.vector(coef(m4, component="stage1")), as.vector(coef(ml4)))
})

test_that("stage1 rank-deficient coefficient covariances are computed correctly with one endog x", {
  expect_equal(vcov(m4, component="stage1"), vcov(ml4))
})

test_that("stage1 rank-deficient coefficients are computed correctly with one endog x, complete=FALSE", {
  expect_equal(as.vector(coef(m4, component="stage1", complete=FALSE)), 
               as.vector(coef(ml4, complete=FALSE)))
})

test_that("stage1 rank-deficient coefficient covariances are computed correctly with one endog x, rank=FALSE", {
  expect_equal(vcov(m4, component="stage1", complete=FALSE), 
               vcov(ml4, complete=FALSE))
})
