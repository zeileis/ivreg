deq <- ivreg(Q ~ P + D | D + F + A, data=Kmenta)
new <- expand.grid(P = seq(90, 110, by=20), D = seq(75, 125, by=20), 
                   F=c(70, 110), A=c(1, 20))
pr <- predict(deq, newdata=new, se.fit=TRUE, interval="prediction")
X <- as.matrix(cbind(1, new[, c("P", "D")]))
yhat <- X %*% coef(deq)
se <- sqrt(diag(X %*% vcov(deq) %*% t(X)))
t <- qt(.025, Inf, lower.tail=FALSE)
lower <- yhat - t*sqrt(se^2 + sigma(deq)^2)
upper <- yhat + t*sqrt(se^2 + sigma(deq)^2)
test_that("predicted values, their std. errors, and prediction intervals are computed correctly", {
  expect_equal(as.vector(cbind(yhat, se, lower, upper)), as.vector(pr))
})
