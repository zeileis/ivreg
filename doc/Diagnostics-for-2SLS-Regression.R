## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, comment=NA, prompt=TRUE, fig.height=4, fig.width=4)

## ------------------------------------------------------------------------
library(ivreg)
Kmenta

## ------------------------------------------------------------------------
deq <- ivreg(Q ~ P + D, ~ D + F + A, data=Kmenta)     # demand equation
summary(deq)

## ------------------------------------------------------------------------
seq <- ivreg(Q ~ P + F + A, ~ D + F + A, data=Kmenta) # supply equation
summary(seq)

## ----fig.height=8, fig.width=8-------------------------------------------
par(mfrow=c(2, 2))
plot(deq)

## ------------------------------------------------------------------------
library(car) # for diagnostic generic functions
qqPlot(deq)

## ------------------------------------------------------------------------
influencePlot(deq)

## ------------------------------------------------------------------------
Kmenta1 <- Kmenta
Kmenta1[20, "Q"] <- 95

## ------------------------------------------------------------------------
deq1 <- update(deq, data=Kmenta1)
compareCoefs(deq, deq1)

## ------------------------------------------------------------------------
qqPlot(deq1)

## ------------------------------------------------------------------------
outlierTest(deq1)

## ------------------------------------------------------------------------
influencePlot(deq1)

## ----fig.height=4, fig.width=8-------------------------------------------
avPlots(deq1)

## ------------------------------------------------------------------------
deq1.20 <- update(deq1, subset = -20)
compareCoefs(deq, deq1, deq1.20)

## ----fig.height=6, fig.width=6-------------------------------------------
H <- cbind(hatvalues(deq1), hatvalues(deq1, type="both"), 
           hatvalues(deq1, type="maximum"))
colnames(H) <- c("stage2", "geom.mean", "maximum")
head(H)
scatterplotMatrix(H, smooth=FALSE)

## ------------------------------------------------------------------------
cbind(dfbeta(deq1)[20, ], coef(deq1) - coef(deq1.20))

## ------------------------------------------------------------------------
c(influence(deq1)$sigma[20], sigma(deq1.20))

## ----fig.height=4, fig.width=8-------------------------------------------
crPlots(deq, smooth=list(span=1))

## ----fig.height=4, fig.width=8-------------------------------------------
library(effects)
plot(predictorEffects(deq, residuals=TRUE), 
     partial.residuals=list(span=1))

## ----fig.height=4, fig.width=8-------------------------------------------
deq2 <- update(deq, . ~ I((P - 85)^4/10^5) + D)
crPlots(deq2, smooth=list(span=1))

## ----fig.height=4, fig.width=8-------------------------------------------
plot(predictorEffects(deq2, residuals=TRUE), 
     partial.residuals=list(span=1))

## ------------------------------------------------------------------------
plot(fitted(deq), rstudent(deq))
abline(h=0)

## ------------------------------------------------------------------------
spreadLevelPlot(deq, smooth=list(span=1))

## ------------------------------------------------------------------------
with(Kmenta, plot(Q, Q^2.5))
abline(lm(Q^2.5 ~ Q, data=Kmenta))

## ------------------------------------------------------------------------
ncvTest(deq)
ncvTest(deq, var = ~ P + D)

## ------------------------------------------------------------------------
summary(deq, vcov=sandwich::sandwich)
SEs <- round(cbind(sqrt(diag(sandwich::sandwich(deq))), 
                   sqrt(diag(vcov(deq)))), 
             4)
colnames(SEs) <- c("sandwich", "conventional")
SEs

## ------------------------------------------------------------------------
Kmenta2 <- Kmenta[, c("D", "F", "A")]
set.seed(492365) # for reproducibility
Kmenta2 <- within(Kmenta2, {
    EQ <- 75.25 + 0.1125*D + 0.1250*F + 0.225*A
    EP <- 85.00 + 0.7500*D - 0.5000*F - 0.900*A
    d1 <- rnorm(20)
    d2 <- rnorm(20)
    v1 <- 2*d1
    v2 <- -0.5*v1 + d2
    v1 <- v1*3*(EQ - min(EQ) + 0.1)/(max(EQ) - min(EQ)) 
              # inducing nonconstant variance
    Q <- EQ + v1
    P <- EP + v2
})

## ------------------------------------------------------------------------
with(Kmenta2, plot(EQ, v1))

## ------------------------------------------------------------------------
deq2 <- update(deq, data=Kmenta2)
summary(deq2)

## ------------------------------------------------------------------------
spreadLevelPlot(deq2)

## ------------------------------------------------------------------------
ncvTest(deq2)

## ------------------------------------------------------------------------
sqrt(vif(deq))

