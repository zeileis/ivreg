## packages
library("ivreg")
library("lmtest")
library("sandwich")

## data and models: full vs. trivial
data("CigaretteDemand", package = "ivreg")
m <- ivreg(log(packs) ~ log(rincome) | log(rprice) | salestax, data = CigaretteDemand)
m0 <- ivreg(log(packs) ~ 1, data = CigaretteDemand)

## heteroscedasticity-consistent HC1 covariance _matrix_ on full model
vc1 <- vcovHC(m, type = "HC1")

## custom HC1 covariance _function_
hc1 <- function(object, ...) vcovHC(object, type = "HC1", ...)

## confint
test_that("confidence intervals are equal with vcov function and matrix", {
  expect_equal(confint(m, vcov = hc1), confint(m, vcov = vc1))
})
test_that("confidence intervals are equal with vcov function and function + ...", {
  expect_equal(confint(m, vcov = hc1), confint(m, vcov = vcovHC, type = "HC1"))
})

## coeftest
test_that("confidence intervals are equal with vcov function and matrix", {
  expect_equal(coeftest(m, vcov = hc1), coeftest(m, vcov = vc1))
})
test_that("coefficient tests are equal with vcov function and function + ...", {
  expect_equal(coeftest(m, vcov = hc1), coeftest(m, vcov = vcovHC, type = "HC1"))
})

## anova
test_that("analysis of variance is equal with vcov function and matrix", {
  expect_equal(anova(m0, m, vcov = hc1), anova(m0, m, vcov = vc1))
})
test_that("analysis of variance is equal with vcov function and function + ...", {
  expect_equal(anova(m0, m, vcov = hc1), anova(m0, m, vcov = vcovHC, type = "HC1"))
})

## summary
test_that("model summary (w/o diagnostics) is equal with vcov function and matrix", {
  expect_equal(summary(m, vcov = hc1, diagnostics = FALSE), summary(m, vcov = vc1, diagnostics = FALSE))
})
test_that("model summary (w/o diagnostics) is equal with vcov function and function + ...", {
  expect_equal(summary(m, vcov = hc1, diagnostics = FALSE), summary(m, vcov = vcovHC, type = "HC1", diagnostics = FALSE))
})
test_that("model summary is equal with vcov function and function + ...", {
  expect_equal(summary(m, vcov = hc1), summary(m, vcov = vcovHC, type = "HC1"))
})
