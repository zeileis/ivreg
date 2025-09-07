m <- ivreg(Q ~ P + D | D + F + A, data=Kmenta)
m1 <- update(m, subset = -1)
test_that("dfbeta computed correctly 1", {
  expect_equal(dfbeta(m)[1, ], coef(m) - coef(m1))
})

m.ivreg <- expect_warning(ivreg(Q ~ P + D | P + D, data=Kmenta)) # OLS
m.lm <- lm(Q ~ P + D, data=Kmenta)

test_that("hatvalues computed correctly", {
  expect_equal(as.vector(hatvalues(m.ivreg)), as.vector(hatvalues(m.lm)))
})

test_that("Cook's distances computed correctly", {
  expect_equal(cooks.distance(m.ivreg), cooks.distance(m.lm))
})

test_that("dfbeta computed correctly 2", {
  expect_equal(dfbeta(m.ivreg), dfbeta(m.lm))
})

Kmenta2 <- Kmenta
Kmenta2$Q[c(1, 4, 10)] <- NA
m.miss <- update(m, data=Kmenta2, na.action=na.exclude)

test_that("deletion statistics computed correctly with na.exclude", {
  expect_true(all(is.na(hatvalues(m.miss)[c(1, 4, 10)])))
  expect_true(all(is.na(residuals(m.miss)[c(1, 4, 10)])))
  expect_true(all(is.na(fitted(m.miss)[c(1, 4, 10)])))
  expect_true(all(is.na(rstudent(m.miss)[c(1, 4, 10)])))
  expect_true(all(is.na(dffits(m.miss)[c(1, 4, 10)])))
  expect_true(all(is.na(cooks.distance(m.miss)[c(1, 4, 10)])))
  expect_true(all(is.na(dfbeta(m.miss)[c(1, 4, 10), ])))
})

m.miss.2<- update(m, data=Kmenta2)
nms <- rownames(Kmenta)[-c(1, 4, 10)]
test_that("rownames of deletion statistics preserved with na.omit", {
  expect_equal(names(hatvalues(m.miss.2)), nms)
  expect_equal(names(residuals(m.miss.2)), nms)
  expect_equal(names(fitted(m.miss.2)), nms)
  expect_equal(names(rstudent(m.miss.2)), nms)
  expect_equal(names(dffits(m.miss.2)), nms)
  expect_equal(names(cooks.distance(m.miss.2)), nms)
  expect_equal(rownames(dfbeta(m.miss.2)), nms)
})

m.ivreg.w <- expect_warning(ivreg(Q ~ P + D | P + D, weights=Q, data=Kmenta)) # WLS
m.lm.w <- lm(Q ~ P + D, data=Kmenta, weights=Q)

test_that("hatvalues computed correctly with weights", {
  expect_equal(as.vector(hatvalues(m.ivreg.w)), as.vector(hatvalues(m.lm.w)))
})

test_that("Cook's distances computed correctly with weights", {
  expect_equal(cooks.distance(m.ivreg.w), cooks.distance(m.lm.w))
})

test_that("rstudent computed correctly", {
  expect_equal(rstudent(m.ivreg), rstudent(m.lm))
})

mw <- ivreg(Q ~ P + D | D + F + A, weights=Q, data=Kmenta)
m1w <- update(mw, subset = -10)
test_that("dfbeta computed correctly with weights", {
  expect_equal(dfbeta(mw)[10, ], coef(mw) - coef(m1w))
})

test_that("rstudent computed correctly with weights", {
  expect_equal(rstudent(m.ivreg.w), rstudent(m.lm.w))
})

test_that("influence measures computed correctly in paralllel", {
  expect_equal(influence(m), influence(m, ncores=2))
})

