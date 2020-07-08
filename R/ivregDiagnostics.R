na.remove <- function(x){
  # remove NAs preserving names
  x <- x[!is.na(x)]
}

diagprod <- function(d, X){
  # equivalent to diag(d) %*% X
  if (!is.vector(d)) stop("d is not a vector")
  if (!is.matrix(X)) stop("X is not a matrix")
  if (length(d) != nrow(X)) stop("d and X not conformable")
  d*X
}

# formula.ivreg <- function(x, ...) formula(x$terms$regressors)

#' Deletion and Other Diagnostic Methods for \code{"ivreg"} Objects
#'
#' @aliases ivregDiagnostics influence.ivreg rstudent.ivreg cooks.distance.ivreg 
#' dfbeta.influence.ivreg hatvalues.ivreg rstudent.influence.ivreg hatvalues.influence.ivreg
#' cooks.distance.influence.ivreg qqPlot.ivreg qqPlot.influence.ivreg influencePlot.ivreg
#' nfluencePlot.influence.ivreg infIndexPlot.ivreg infIndexPlot.influence.ivreg
#' model.matrix.influence.ivreg avPlot.ivreg Boot.ivreg 
#' crPlot.ivreg plot.ivreg qqPlot.ivreg outlierTest.ivreg influencePlot.ivreg 
#' spreadLevelPlot.ivreg ncvTest.ivreg deviance.ivreg
#' 
#' @param model,x,object A \code{"ivreg"} or \code{"influence.ivreg"} object.
#' @param sigma. If \code{TRUE} (the default for 1000 or fewer cases), the deleted value
#' of the residual standard deviation is computed for each case; if \code{FALSE}, the
#' overall residual standard deviation is used to compute other deletion diagnostics.
#' @param applyfun Optional loop replacement function that should work like
#' \code{\link[base]{lapply}} with arguments \code{function(X, FUN, ...)}. The default
#' is to use a loop unless the \code{ncores} argument is specified (see below).
#' @param ncores Numeric, number of cores to be used in parallel computations. If set
#' to an integer the \code{applyfun} is set to use either \code{\link[parallel]{parLapply}}
#' (on Windows) or \code{\link[parallel]{mclapply}} (otherwise) with the desired number of cores.
#' @param type If \code{"stage2"} (the default), hatvalues are for the second stage regression;
#' if \code{"both"}, the hatvalues are the geometric mean of the casewise hatvalues for the
#' two stages; if \code{"maximum"}, the hatvalues are the larger of the casewise
#' hatvalues for the two stages. In computing the geometric mean or casewise maximum hatvalues,
#' the hatvalues for each stage are first divided by their average (number of coefficients in
#' stage regression/number of cases); the geometric mean or casewise maximum values are then
#' multiplied by the average hatvalue from the second stage.
#' @param main Main title for the graph.
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
#' @description Methods for computing deletion and other regression diagnostics for 2SLS regression.
#' It's generally more efficient to compute the deletion diagnostics via the \code{influence}
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
#' @seealso \code{\link{ivreg}}, \code{\link[car]{avPlots}},
#'   \code{\link[car]{crPlots}}, \code{\link[effects]{predictorEffects}},
#'   \code{\link[car]{qqPlot}}, \code{\link[car]{influencePlot}},
#'   \code{\link[car]{infIndexPlot}}, \code{\link[car]{Boot}}, 
#'   \code{\link[car]{outlierTest}}, \code{\link[car]{spreadLevelPlot}}, 
#'   \code{\link[car]{ncvTest}}.
#' @examples
#' kmenta.eq1 <- ivreg(Q ~ P + D | D + F + A, data = Kmenta)
#' summary(kmenta.eq1)
#' car::avPlots(kmenta.eq1)
#' car::crPlots(kmenta.eq1)
#' car::influencePlot(kmenta.eq1)
#' car::influenceIndexPlot(kmenta.eq1)
#' car::qqPlot(kmenta.eq1)
#' car::spreadLevelPlot(kmenta.eq1)
#' plot(effects::predictorEffects(kmenta.eq1, residuals = TRUE))
#' set.seed <- 12321 # for reproducibility
#' confint(car::Boot(kmenta.eq1, R = 250)) # 250 reps for brevity
#' car::outlierTest(kmenta.eq1)
#' car::ncvTest(kmenta.eq1)
#' 
#' @rdname ivregDiagnostics
#' @importFrom stats influence
#' @export
influence.ivreg <- function(model, sigma. = n <= 1e3, type = c("stage2", "both", "maximum"), 
                            applyfun = NULL, ncores = NULL, ...){

  type <- match.arg(type)

  Z <- model.matrix(model, component = "instruments") # model$model.matrix.instruments
  X <- model.matrix(model, component = "regressors") # model$model.matrix
  X.fit <- model.matrix(model, component = "projected") # model$fitted1 
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
  
  ## set up parallel apply if specified
  if(!is.null(ncores) && is.null(applyfun)) {
    applyfun <- if(ncores == 1L) {
      lapply
    } else if(.Platform$OS.type == "windows") {
      cl <- parallel::makeCluster(ncores) 
      on.exit(parallel::stopCluster(cl))
      function(X, FUN, ...) parallel::parLapply(cl, X, FUN, ...)
    } else {
      function(X, FUN, ...) parallel::mclapply(X, FUN, ..., mc.cores = ncores)    
    }
  }

  if (is.function(applyfun)){
    result <- applyfun(1:n, function(i) {
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
    })
    result <- do.call("rbind", result)
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
                 coefficients = coef(model),
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

#' @rdname ivregDiagnostics
#' @importFrom stats rstudent
#' @export
rstudent.ivreg <- function(model, ...) {
  influence(model)$rstudent
}

#' @rdname ivregDiagnostics
#' @importFrom stats cooks.distance
#' @method cooks.distance ivreg
#' @export
cooks.distance.ivreg <- function(model, ...) {
  influence(model)$cookd
}

#' @rdname ivregDiagnostics
#' @export
dfbeta.influence.ivreg <- function(model, ...) {
  model$dfbeta
}

#' @rdname ivregDiagnostics
#' @importFrom stats dfbeta
#' @export
dfbeta.ivreg <- function(model, ...) {
  influence(model)$dfbeta
}

#' @rdname ivregDiagnostics
#' @importFrom stats df.residual hatvalues lm.influence na.omit naresid
#' @export
hatvalues.ivreg <- function(model, type = c("stage2", "both", "maximum"), ...){
  type <- match.arg(type)
  # if (!inherits(model, "lm")) {
  #   class(model) <- c(class(model), "lm")
  #   return(hatvalues(model, type=type, ...))
  # }
  hatvalues <- if (type == "stage2") {
    .Class <- "lm"
    NextMethod()
  } else {
    n <- model$nobs
    p <- model$p
    q <- model$q
    mean1 <- q/n
    mean2 <- p/n
    hat2 <- lm.influence(model)$hat/mean2
    model[c("qr", "rank", "residuals", "coefficients")] <-
      list(model$qr1, model$rank1, na.omit(residuals(model, type="projected")), 
           model$coefficients1)
    hat1 <- lm.influence(model)$hat/mean1
    hat <- if (type == "both") {
      sqrt(hat1*hat2)
    } else {
      pmax(hat1, hat2)
    }
    mean2*hat
  }
  na.action <- model$na.action
  if(inherits(na.action, "exclude")) hatvalues[na.action]  <- NA
  hatvalues
}

#' @rdname ivregDiagnostics
#' @method rstudent influence.ivreg
#' @export
rstudent.influence.ivreg <- function(model, ...) {
  model$rstudent
}

#' @rdname ivregDiagnostics
#' @method hatvalues influence.ivreg
#' @export
hatvalues.influence.ivreg <- function(model, ...) {
  model$hatvalues
}

#' @rdname ivregDiagnostics
#' @export
#' @method cooks.distance influence.ivreg
cooks.distance.influence.ivreg <- {
  function(model, ...) model$cookd
}

##' @importFrom car qqPlot
##' @importFrom graphics par
##' @export
# qqPlot.ivreg <- function(x,
#                          ylab = paste("Studentized Residuals(", deparse(substitute(x)), ")", sep = ""),
#                          distribution = c("t", "norm"), ...){
#   distribution <- match.arg(distribution)
#   rstudent <- rstudent(x)
#   if (distribution == "t"){
#     car::qqPlot(rstudent, ylab = ylab, distribution = "t", df = df.residual(x), ...)
#   } else {
#     car::qqPlot(rstudent, ylab = ylab, distribution = "norm", ...)
#   }
# }

#' @rdname ivregDiagnostics
#' @method qqPlot influence.ivreg
#' @param distribution \code{"t"} (the default) or \code{"norm"}.
#' @param ylab The vertical axis label.
#' @export
qqPlot.influence.ivreg <- function(x,
                                   ylab = paste("Studentized Residuals(", deparse(substitute(x)), ")", sep = ""),
                                   distribution = c("t", "norm"), ...){
  distribution <- match.arg(distribution)
  rstudent <- rstudent(x)
  if (distribution == "t"){
    car::qqPlot(rstudent, ylab = ylab, distribution = "t", df = df.residual(x), ...)
  } else {
    car::qqPlot(rstudent, ylab = ylab, ...)
  }
}

#' @rdname ivregDiagnostics
#' @method influencePlot ivreg
#' @importFrom car influencePlot
#' @export
influencePlot.ivreg <- function(model, ...){
  influencePlot(influence(model), ...)
}

#' @rdname ivregDiagnostics
#' @method influencePlot influence.ivreg
#' @export
influencePlot.influence.ivreg <- function(model, ...){
  # if (!inherits(model, "lm")) {
  #   class(model) <- c(class(model), "lm")
  #   influencePlot(model)
  # }
  # else NextMethod()
  .Class <- "lm"
  NextMethod()
}

#' @rdname ivregDiagnostics
#' @method infIndexPlot ivreg
#' @export
infIndexPlot.ivreg <- function(model, ...){
  infIndexPlot(influence(model), ...)
}

#' @rdname ivregDiagnostics
#' @method infIndexPlot influence.ivreg
#' @importFrom car infIndexPlot
#' @export
infIndexPlot.influence.ivreg <- function(model, ...){
  if (!inherits(model, "lm")) {
    class(model) <- c(class(model), "lm")
    infIndexPlot(model, ...)
  }
  else NextMethod()
  # .Class <- "lm"
  # NextMethod()
}

#' @rdname ivregDiagnostics
#' @method model.matrix influence.ivreg
#' @export
model.matrix.influence.ivreg <- function(object, ...){
  object$model
}

#' @rdname ivregDiagnostics
#' @importFrom car avPlot
#' @export
avPlot.ivreg <- function(model, ...){
  # if (!inherits(model, "lm")) {
  #   class(model) <- c(class(model), "lm")
  #   model$model.matrix <- model.matrix(model, type = "projected")
  #   avPlot(model, ...)
  # } else {
  #   NextMethod()
  # }
  .Class <- "lm"
  NextMethod()
}

#' @rdname ivregDiagnostics
#' @method Boot ivreg
#' @importFrom car Boot
#' @importFrom stats coef
#' @param method only \code{"case"} (case resampling) is supported: see \code{\link[car]{Boot}}.
#' @param f,labels,R see \code{\link[car]{Boot}}.
#' @export
Boot.ivreg <- function(object, f = coef, labels = names(f(object)), R = 999, 
                       method = "case", ncores = 1, ...){
  method <- match.arg(method)
  NextMethod()
}

#' @rdname ivregDiagnostics
#' @importFrom car crPlot
#' @export
crPlot.ivreg <- function(model, ...){
  # if (!inherits(model, "lm")) {
  #   class(model) <- c(class(model), "lm")
  #   crPlot(model, ...)
  # } else {
  #   NextMethod()
  # }
  .Class <- "lm"
  NextMethod()
}

#' @rdname ivregDiagnostics
#' @importFrom graphics plot
#' @export
plot.ivreg <- function(x, ...){
  if (!inherits(x, "lm")) {
    class(x) <- c(class(x), "lm")
    plot(x, ...)
  } else {
    NextMethod()
  }
  # .Class <- "lm"
  # NextMethod()
}

#' @rdname ivregDiagnostics
#' @importFrom car qqPlot
#' @export
qqPlot.ivreg <- function(x, distribution=c("t", "norm"), ...){
  distribution <- match.arg(distribution)
  rstudent <- rstudent(x)
  df <- df.residual(x)
  qqPlot(rstudent, distribution=distribution, df=df, ylab="Studentized Residuals", ...)
}

#' @rdname ivregDiagnostics
#' @importFrom car outlierTest
#' @export
outlierTest.ivreg <- function(x, ...){
  # if (!inherits(x, "lm")) {
  #   class(x) <- c(class(x), "lm")
  #   outlierTest(x, ...)
  # } else {
  #   NextMethod()
  # }
  .Class <- "lm"
  NextMethod()
}

#' @rdname ivregDiagnostics
#' @importFrom car influencePlot
#' @export
influencePlot.ivreg <- function(x, ...){
  # if (!inherits(x, "lm")) {
  #     class(x) <- c(class(x), "lm")
  #   influencePlot(x, ...)
  # } else {
  #   NextMethod()
  # }
  .Class <- "lm"
  NextMethod()
}

#' @rdname ivregDiagnostics
#' @importFrom car spreadLevelPlot
#' @export
spreadLevelPlot.ivreg <- function(x, main="Spread-Level Plot", ...){
  # if (!inherits(x, "lm")) {
  #   class(x) <- c(class(x), "lm")
  #   spreadLevelPlot(x, main=main, ...)
  # } else {
  #   NextMethod()
  # }
  .Class <- "lm"
  NextMethod()
}

#' @rdname ivregDiagnostics
#' @importFrom car ncvTest
#' @export
ncvTest.ivreg <- function(model, ...){
  # if (!inherits(model, "lm")) {
  #   class(model) <- c(class(model), "lm")
  #   ncvTest(model, ...)
  # } else {
  #   NextMethod()
  # }
    .Class <- "lm"
    NextMethod()
}

#' @rdname ivregDiagnostics
#' @importFrom stats deviance
#' @export
deviance.ivreg <- function(object, ...){
  # if (!inherits(object, "lm")) {
  #   class(object) <- c(class(object), "lm")
  #   deviance(object, ...)
  # } else {
  #   NextMethod()
  # }
  .Class <- "lm"
  NextMethod()
}
