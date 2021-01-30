if (require(MASS)){
  
  Kmenta1 <- Kmenta
  Kmenta1[20, "Q"] <- 95 # corrupted data
  deqr1 <- ivreg(Q ~ P + D | D + F + A, data=Kmenta1, method="MM") 
  stg1 <- rlm(P ~ D + F + A, data=Kmenta1, method="MM")
  Kmenta1$P.hat <- fitted(stg1)
  stg2 <- rlm(Q ~ P.hat + D, data=Kmenta1, method="MM")
  
  test_that("robust 2SLS coefficients are computed correctly", {
    expect_equal(coef(deqr1), coef(stg2), check.attributes=FALSE)
  })
  
  test_that("hatvalues are computed correctly for robust 2SLS, stage 2",{
    hats <- c(hatvalues(lm(residuals(stg2) ~ model.matrix(stg2) - 1, weights=stg2$w)), 0)
    names(hats)[20] <- "1941"
    expect_equal(hats, hatvalues(deqr1))
  })
  
  test_that("hatvalues are computed correctly for robust 2SLS, stage 1",{
    expect_equal(hatvalues(stg1), hatvalues(deqr1, type="stage1"))
  })
  
}
