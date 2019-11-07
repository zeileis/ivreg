na.remove <- function(x){
  # remove NAs preserving names
  x <- x[!is.na(x)]
}

formula.ivreg <- function(x, ...) formula(x$terms$regressors)

#' Deletion and Other Diagnostic Methods for \code{"ivreg"} Objects
#'
#' @aliases 2SLS_Diagnostics influence.ivreg rstudent.ivreg cooks.distance.ivreg 
#' dfbeta.influence.ivreg hatvalues.ivreg rstudent.influence.ivreg hatvalues.influence.ivreg
#' cooks.distance.influence.ivreg qqPlot.ivreg qPlot.influence.ivreg influencePlot.ivreg
#' nfluencePlot.influence.ivreg infIndexPlot.ivreg infIndexPlot.influence.ivreg
#' model.matrix.influence.ivreg avPlot.ivreg Boot.ivreg
#' 
#' @param model,x,object A \code{"ivreg"} or \code{"influence.ivreg"} object.
#' @param sigma. If \code{TRUE} (the default for 1000 or fewer cases), the deleted value
#' of the residual standard deviation is computed for each case; if \code{FALSE}, the
#' overall residual standard deviation is used to compute other deletion diagnostics.
#' @param ncores If \code{ncores > 1} (the default is \code{1}) computations are performed
#' in parallel; this is advantageous only for large data sets.
#' @param type If \code{"stage2"} (the default), hatvalues are for the second stage regression;
#' if \code{"both"}, the hatvalues are the geometric mean of the casewise hatvalues for the
#' two stages; if \code{"maximum"}, the hatvalues are the larger of the casewise
#' hatvalues for the two stages. In computing the geometric mean or casewise maximum hatvalues,
#' the hatvalues for each stage are first divided by their average (number of coefficients in
#' stage regression/number of cases); the geometric mean or casewise maximum values are then
#' multiplied by the average hatvalue from the second stage.
#' @param ... arguments to be passed down.
#'
#' @return In the case of \code{influence.ivreg}, an object of class \code{"influence.ivreg"}
#' with the following components:
#' \describe{
#' \item{\code{coefficients}}{the estimated regression coefficients}
#' \item{\code{model}}{the model matrix}
#' \item{\code{dfbeta}}{influence on coefficients}
#' \item{\code{sigma}}{deleted values of the residual standard deviation}
#' \item{\code{dffits}}{overall influence on the regression coefficients}
#' \item{\code{cookd}}{Cook's distances}
#' \item{\code{hatvalues}}{hatvalues}
#' \item{\code{rstudent}}{Studentized residuals}
#' \item{\code{df.residual}}{residual degrees of freedom}
#' }
#' In the case of other methods, such as \code{rstudent.ivreg} or
#' \code{rstudent.influence.ivreg}, the corresponding diagnostic statistics.
#'
#' @description Methods for computing deletion diagnostics for 2SLS regression.
#' It's generally more efficient to compute the diagnostics via the \code{influence}
#' method and then to extract the various specific diagnostics with the methods for
#' \code{"influence.ivreg"} objects. Other diagnostics for linear models, such as
#' added-variable plots (\code{\link[car]{avPlots}}) and component-plus-residual
#' plots (\code{\link[car]{crPlots}}), also work, as do effect plots
#' (e.g., \code{\link[effects]{predictorEffects}}) with residuals (see the examples below).
#' The pointwise confidence envelope for the \code{\link[car]{qqPlot}} method assumes an independent random sample
#' from the t distribution with degrees of freedom equal to the residual degrees of
#' freedom for the model and so are approximate, because the studentized residuals aren't
#' independent.
#' 
#' For additional information, see the vignette 
#' \href{../doc/Diagnostics-for-2SLS-Regression.pdf}{Diagnostics for 2SLS Regression}.
#' @import foreach
#' @importFrom stats influence
#' @export
#' @seealso \code{\link{ivreg}}, \link{2SLS_Methods}, \code{\link[car]{avPlots}},
#'   \code{\link[car]{crPlots}}, \code{\link[effects]{predictorEffects}},
#'   \code{\link[car]{qqPlot}}, \code{\link[car]{influencePlot}},
#'   \code{\link[car]{infIndexPlot}}, \code{\link[car]{Boot}}.
#' @examples
#' kmenta.eq1 <- ivreg(Q ~ P + D | D + F + A, data=Kmenta)
#' car::avPlots(kmenta.eq1)
#' car::crPlots(kmenta.eq1)
#' car::influencePlot(kmenta.eq1)
#' car::influenceIndexPlot(kmenta.eq1)
#' car::qqPlot(kmenta.eq1)
#' plot(effects::predictorEffects(kmenta.eq1, residuals=TRUE))
#' set.seed <- 12321 # for reproducibility
#' confint(car::Boot(kmenta.eq1, R=250)) # 250 reps for brevity
influence.ivreg <- function(model, sigma. = n <= 1e3, type=c("stage2", "both", "maximum"), 
                            ncores=1, ...){

  type <- match.arg(type)

  Z <- model.matrix(model, component="instruments") # model$model.matrix.instruments
  X <- model.matrix(model, component="regressors") # model$model.matrix
  X.fit <- model.matrix(model, component="projected") # model$fitted.1 
  y <- model$y
  if (is.null(y)) stop("response variable not in model object")
  b <- coef(model) # model$coefficients
  res <- na.remove(residuals(model)) # na.remove(model$residuals)
  .sigma <- sqrt(model$sigma^2)
  hatvalues <-  na.remove(hatvalues(model, type=type))

  na.action <- model$na.action 

  rnames <- rownames(X)
  cnames <- colnames(X)

  names(hatvalues) <- rnames

  w <- na.remove(weights(model)) # na.remove(model$weights)
  if (!is.null(w)){
    w <- sqrt(w)
    X <- diagprod(w, X)
    Z <- diagprod(w, Z)
    X.fit <- diagprod(w, X.fit)
    y <- w*y
  }
  else w <- 1

  rss <- sum((w*res)^2)
  ZtZinv <- solve(crossprod(Z)) #TODO: avoid matrix inversions?
  XtZ <- crossprod(X, Z)
  A <- XtZ %*% ZtZinv %*% t(XtZ)
  Ainv <- solve(A)
  pi <- ZtZinv %*% crossprod(Z, y)
  r <- XtZ %*% ZtZinv %*% t(Z)
  XfXfinv <- solve(crossprod(X.fit))

  n <- model$nobs
  p <- length(b)  # model$p
  
  if (ncores > 1){
    cl <- parallel::makeCluster(ncores) 
    doParallel::registerDoParallel(cl)
    result <- foreach(i = 1:n, .combine=rbind) %dopar% {
      c <- as.vector(Z[i, ] %*% ZtZinv %*% Z[i, ])
      Xmr <- X[i, ] - r[, i]
      XiAinvXi <- as.vector(X[i, ] %*% Ainv %*% X[i, ])
      XmrAinvXi <- as.vector(Xmr %*% Ainv %*% X[i, ])
      XmrAinvXmr <- as.vector(Xmr %*% Ainv %*% Xmr)
      delta <- 1 - XiAinvXi + XmrAinvXi^2/(1 - c + XmrAinvXmr)
      denom <- (1 - c + XmrAinvXmr)*delta
      h <- Xmr * (1 - XiAinvXi) / denom + X[i, ] * XmrAinvXi / denom
      j <- Xmr * XmrAinvXi / denom - X[i, ]/delta
      g <- h * as.vector((y[i] - Z[i, ] %*% pi) - (y[i] - r[, i] %*% b)) +
        (h + j)*as.vector(y[i] - X[i, ] %*% b)
      dfbeta.i <- - as.vector(Ainv %*% g)
      sigma.i <- if (sigma.){
        ss <- rss + as.vector(g %*% Ainv %*% crossprod(X[-i, ]) %*% Ainv %*% g) -
          2 * as.vector(g %*% Ainv %*% t(X[-i, ]) %*% (y[-i] - X[-i, ] %*% b))  -
          as.vector(y[i] - X[i, ] %*% b)^2
        sqrt(ss/(n - p - 1))
      } else .sigma
      dffits.i <- as.vector(X[i, ] %*% dfbeta.i)/
        (sigma.i * as.vector(sqrt(X[i, ] %*% XfXfinv %*% X[i, ])))
      c(dfbeta.i, sigma.i, dffits.i)
    }
    parallel::stopCluster(cl)
    nc <- ncol(result)
    dfbeta <- result[, -c(nc - 1, nc)]
    sigma <- result[, nc - 1]
    dffits <- result[, nc]
  } else {
    dfbeta <- matrix(0, n, p)
    dffits <- cookd <- rep(0, n)
    sigma <- rep(.sigma, n)
    for (i in 1:n){ #TODO: move this loop to cpp code?
      c <- as.vector(Z[i, ] %*% ZtZinv %*% Z[i, ])
      Xmr <- X[i, ] - r[, i]
      XiAinvXi <- as.vector(X[i, ] %*% Ainv %*% X[i, ])
      XmrAinvXi <- as.vector(Xmr %*% Ainv %*% X[i, ])
      XmrAinvXmr <- as.vector(Xmr %*% Ainv %*% Xmr)
      delta <- 1 - XiAinvXi + XmrAinvXi^2/(1 - c + XmrAinvXmr)
      denom <- (1 - c + XmrAinvXmr)*delta
      h <- Xmr * (1 - XiAinvXi) / denom + X[i, ] * XmrAinvXi / denom
      j <- Xmr * XmrAinvXi / denom - X[i, ]/delta
      g <- h * as.vector((y[i] - Z[i, ] %*% pi) - (y[i] - r[, i] %*% b)) +
        (h + j)*as.vector(y[i] - X[i, ] %*% b)
      dfbeta[i, ] <- - Ainv %*% g
      if (sigma.){
        ss <- rss + as.vector(g %*% Ainv %*% crossprod(X[-i, ]) %*% Ainv %*% g) -
          2 * as.vector(g %*% Ainv %*% t(X[-i, ]) %*% (y[-i] - X[-i, ] %*% b))  -
          as.vector(y[i] - X[i, ] %*% b)^2
        sigma[i] <- sqrt(ss/(n - p - 1))
      }
      dffits[i] <- as.vector(X[i, ] %*% dfbeta[i, ])/
        (sigma[i] * as.vector(sqrt(X[i, ] %*% XfXfinv %*% X[i, ])))
    }
  }

  rstudent <- w*res/(sigma * sqrt(1 - hatvalues))
  cookd <- (sigma^2/.sigma^2)*dffits^2/p
  
  rownames(dfbeta) <- rnames
  colnames(dfbeta) <- cnames
  names(dffits) <- names(sigma) <- names(cookd) <- rnames

  result <- list(model = model.matrix(model),
                 coefficients=coef(model),
                 dfbeta = naresid(na.action, dfbeta),
                 sigma = naresid(na.action, sigma),
                 dffits = naresid(na.action, dffits),
                 cookd = naresid(na.action, cookd),
                 hatvalues = naresid(na.action, hatvalues),
                 rstudent = naresid(na.action, rstudent),
                 df.residual = df.residual(model))
  class(result) <- "influence.ivreg"
  result
}

#' @rdname influence.ivreg
#' @importFrom stats rstudent
#' @export
rstudent.ivreg <- function(model, ...) {
  influence(model)$rstudent
}

#' @rdname influence.ivreg
#' @importFrom stats cooks.distance
#' @method cooks.distance ivreg
#' @export
cooks.distance.ivreg <- function(model, ...) {
  influence(model)$cookd
}

#' @rdname influence.ivreg
#' @export
dfbeta.influence.ivreg <- function(model, ...) {
  model$dfbeta
}

#' @rdname influence.ivreg
#' @importFrom stats dfbeta
#' @export
dfbeta.ivreg <- function(model, ...) {
  influence(model)$dfbeta
}

#' @rdname influence.ivreg
#' @importFrom stats hatvalues lm.influence
#' @export
hatvalues.ivreg <- function(model, type=c("stage2", "both", "maximum"), ...){
  type <- match.arg(type)
  hatvalues <- if (type == "stage2") NextMethod() else {
    n <- model$nobs
    p <- model$p
    q <- model$q
    mean1 <- q/n
    mean2 <- p/n
    hat.2 <- lm.influence(model)$hat/mean2
    model[c("qr", "rank", "residuals", "coefficients")] <-
      list(model$qr.1, model$rank.1, na.omit(residuals(model, type="projected")), 
           model$coefficients.1)
    hat.1 <- lm.influence(model)$hat/mean1
    hat <- if (type == "both") {
      sqrt(hat.1*hat.2)
    } else {
      pmax(hat.1, hat.2)
    }
    mean2*hat
  }
  na.action <- model$na.action
  if(class(na.action) == "exclude") hatvalues[na.action]  <- NA
  hatvalues
}

#' @rdname influence.ivreg
#' @method rstudent influence.ivreg
#' @export
rstudent.influence.ivreg <- function(model, ...) {
  model$rstudent
}

#' @rdname influence.ivreg
#' @method hatvalues influence.ivreg
#' @export
hatvalues.influence.ivreg <- function(model, ...) {
  model$hatvalues
}

#' @rdname influence.ivreg
#' @export
#' @method cooks.distance influence.ivreg
cooks.distance.influence.ivreg <- {
  function(model, ...) model$cookd
}

#' @rdname influence.ivreg
#' @importFrom car qqPlot
#' @importFrom graphics par
#' @export
qqPlot.ivreg <- function(x,
                        ylab=paste("Studentized Residuals(",deparse(substitute(x)), ")", sep=""),
                        distribution=c("t", "norm"), ...){
  distribution <- match.arg(distribution)
  rstudent <- rstudent(x)
  if (distribution == "t"){
    car::qqPlot(rstudent, ylab=ylab, distribution="t", df=df.residual(x), ...)
  } else {
    car::qqPlot(rstudent, ylab=ylab, distribution="norm", ...)
  }
}

#' @rdname influence.ivreg
#' @method qqPlot influence.ivreg
#' @param distribution \code{"t"} (the default) or \code{"norm"}.
#' @param ylab The vertical axis label.
#' @export
qqPlot.influence.ivreg <- function(x,
                                  ylab=paste("Studentized Residuals(",deparse(substitute(x)), ")", sep=""),
                                  distribution=c("t", "norm"), ...){
  distribution <- match.arg(distribution)
  rstudent <- rstudent(x)
  if (distribution == "t"){
    car::qqPlot(rstudent, ylab=ylab, distribution="t", df=df.residual(x), ...)
  } else {
    car::qqPlot(rstudent, ylab=ylab, ...)
  }
}

#' @rdname influence.ivreg
#' @method influencePlot ivreg
#' @importFrom car influencePlot
#' @export
influencePlot.ivreg <- function(model, ...){
  influencePlot(influence(model), ...)
}

#' @rdname influence.ivreg
#' @method influencePlot influence.ivreg
#' @export
influencePlot.influence.ivreg <- function(model, ...){
  if (length(class(model)) == 1) {
    class(model) <- c(class(model), "lm")
    influencePlot(model)
  }
  else NextMethod()
}

#' @rdname influence.ivreg
#' @method infIndexPlot ivreg
#' @export
infIndexPlot.ivreg <- function(model, ...){
  infIndexPlot(influence(model), ...)
}

#' @rdname influence.ivreg
#' @method infIndexPlot influence.ivreg
#' @importFrom car infIndexPlot
#' @export
infIndexPlot.influence.ivreg <- function(model, ...){
  if (length(class(model)) == 1) {
    class(model) <- c(class(model), "lm")
    infIndexPlot(model, ...)
  }
  else NextMethod()
}

#' @rdname influence.ivreg
#' @method model.matrix influence.ivreg
#' @export
model.matrix.influence.ivreg <- function(object, ...){
  object$model
}

#' @rdname influence.ivreg
#' @importFrom car avPlot
#' @export
avPlot.ivreg <- function(model, ...){
  model$model.matrix <- model.matrix(model, type="projected")
  NextMethod()
}

#' @rdname influence.ivreg
#' @method Boot ivreg
#' @importFrom car Boot
#' @param method only \code{"case"} (case resampling) is supported: see \code{\link[car]{Boot}}.
#' @param f,labels,R see \code{\link[car]{Boot}}.
#' @export
Boot.ivreg <- function(object, f = coef, labels = names(f(object)), R = 999, 
                       method = "case", ncores = 1, ...){
  method <- match.arg(method)
  NextMethod()
}
