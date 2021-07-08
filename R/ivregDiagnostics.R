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
#' influencePlot.influence.ivreg infIndexPlot.ivreg infIndexPlot.influence.ivreg
#' model.matrix.influence.ivreg avPlots.ivreg avPlot.ivreg mcPlots.ivreg mcPlot.ivreg Boot.ivreg
#' crPlots.ivreg crPlot.ivreg ceresPlots.ivreg ceresPlot.ivreg plot.ivreg
#' outlierTest.ivreg influencePlot.ivreg
#' spreadLevelPlot.ivreg ncvTest.ivreg deviance.ivreg
#'
#' @param model,x,object A \code{"ivreg"} or \code{"influence.ivreg"} object.
#' @param sigma. If \code{TRUE} (the default for 1000 or fewer cases), the deleted value
#' of the residual standard deviation is computed for each case; if \code{FALSE}, the
#' overall residual standard deviation is used to compute other deletion diagnostics.
#' @param applyfun Optional loop replacement function that should work like
#' \code{\link[base]{lapply}} with arguments \code{function(X, FUN, ...)}. The default
#' is to use a loop unless the \code{ncores} argument is specified (see below).
#' @param ncores
#' Numeric, number of cores to be used in parallel computations. If set
#' to an integer the \code{applyfun} is set to use either \code{\link[parallel:clusterApply]{parLapply}}
#' (on Windows) or
#' #ifdef windows
#' \code{\link[parallel:mcdummies]{mclapply}}
#' #endif
#' #ifdef unix
#' \code{\link[parallel]{mclapply}}
#' #endif
#' (otherwise) with the desired number of cores.
#'
#' @param type If \code{"stage2"} (the default), hatvalues are for the second stage regression;
#' if \code{"both"}, the hatvalues are the geometric mean of the casewise hatvalues for the
#' two stages; if \code{"maximum"}, the hatvalues are the larger of the casewise
#' hatvalues for the two stages. In computing the geometric mean or casewise maximum hatvalues,
#' the hatvalues for each stage are first divided by their average (number of coefficients in
#' stage regression/number of cases); the geometric mean or casewise maximum values are then
#' multiplied by the average hatvalue from the second stage.
#' @param terms Terms for which added-variable plots are to be constructed; the default,
#' if the argument isn't specified, is the \code{"regressors"} component of the model formula.
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
#' Many other methods (e.g., \code{crPlot.ivreg}, \code{avPlot.ivreg}, \code{Effect.ivreg})
#' draw graphs.
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
#' car::mcPlots(kmenta.eq1)
#' car::crPlots(kmenta.eq1)
#' car::ceresPlots(kmenta.eq1)
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

  type <- match.arg(type, c("stage2", "both", "maximum"))

  Z <- model.matrix(model, component = "instruments") # model$model.matrix.instruments
  X <- model.matrix(model, component = "regressors") # model$model.matrix
  X_proj <- model.matrix(model, component = "projected") # model$fitted1
  y <- model$y
  vcov_unscaled <- model$cov.unscaled
  if (is.null(y)) stop("response variable not in model object")
  b <- coef(model) # model$coefficients
  res <- na.remove(residuals(model)) # na.remove(model$residuals)
  # Residuals of first stage and the projection
  X_resid <- if (is.null(model$residuals1)) X - X_proj else model$residuals1
  res_p <- if (is.null(model$residuals2)) y - X_proj %*% b else model$residuals2
  .sigma <- model$sigma
  w <- na.remove(weights(model)) # na.remove(model$weights)
  hatvalues <-  hatvalues(model, type=type)
  if (!is.null(w) && length(hatvalues) != length(w)){
    h <- rep(0, length(w))
    h[w > 0] <- hatvalues
    hatvalues <- h
  }

  na.action <- model$na.action

  rnames <- rownames(X)
  cnames <- colnames(X)

  names(hatvalues) <- rnames
  hatvalues <- na.remove(hatvalues)

  if (!is.null(w)){
    w <- sqrt(w)
    X <- diagprod(w, X)
    Z <- diagprod(w, Z)
    X_proj <- diagprod(w, X_proj)
    X_resid <- diagprod(w, X_resid)
    res_p <- w * res_p
    res <- w * res
    y <- w * y
  } else w <- 1

  rss <- sum(res^2)
  qr_Z <- qr(Z)
  # Residuals of y ~ Z
  res_z <- qr.resid(qr_Z, y)
  qr_A <- qr(crossprod(X, X_proj)) # A = X' Pz X -- we reuse the QR

  n <- model$nobs
  p <- length(b)  # model$p
  idx <- seq_len(n)

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

  Ai_X <- qr.solve(qr_A, t(X))
  Ai_Xr <- qr.solve(qr_A, t(X_resid))
  h_Pz <- rowSums(qr.Q(qr_Z)^2) # Diagonal of the projection matrix P_z
  # diag(X (A)^-1 X')
  h_X <- vapply(idx, function(i) {X[i, ] %*% Ai_X[, i]}, numeric(1L))
  # diag(X_resid (A)^-1 X')
  h_XrX <- vapply(idx, function(i) {X_resid[i, ] %*% Ai_X[, i]}, numeric(1L))
  # diag(X_resid (A)^-1 X_resid)
  h_XrXr <- vapply(idx, function(i) {X_resid[i, ] %*% Ai_Xr[, i]}, numeric(1L))

  denom <- (1 - h_Pz + h_XrXr)
  # If delta is tiny we may have numerical issues
  delta <- 1 - h_X + (h_XrX^2) / denom
  h <- ((1 - h_X) * X_resid + (h_XrX) * X) / (denom * delta)
  j <- ((h_XrX * X_resid) / denom - X) / delta
  # Now we have DFBETA, loosely following Phillips (1977)
  g <- h * as.numeric(res_z - res_p) + (h + j) * res
  dfbeta <- t(-qr.solve(qr_A, t(g)))

  if(isTRUE(sigma.)) {
    dfb_sq <- vapply(idx, function(i) {
      -dfbeta[i, ] %*% crossprod(X[-i, ]) %*% -dfbeta[i, ]}, numeric(1L))
    dfb_resid <- vapply(idx, function(i) {
      -dfbeta[i, ] %*% crossprod(X[-i, ], res[-i])}, numeric(1L))
    rss_i <- rss + dfb_sq - 2 * dfb_resid - res^2
    sigma <- sqrt(rss_i / (n - p - 1L))
  } else {
    sigma <- rep(.sigma, n)
  }
  dffits <- vapply(idx, function(i) {
    X[i, ] %*% dfbeta[i, ] / (sigma[i] * sqrt(h_X[i]))}, numeric(1L))

  rstudent <- res / (sigma * sqrt(1 - naresid(na.action, hatvalues)))
  cookd <- (sigma^2 / .sigma^2) * dffits^2 / p

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
hatvalues.ivreg <- function(model, type = c("stage2", "both", "maximum", "stage1"), ...){
  type <- match.arg(type, c("stage2", "both", "maximum", "stage1"))
  hats <- model$hatvalues
  nms <- names(model$residuals)
  hatvalues <- if (type == "stage2") {
    # .Class <- "lm"
    # NextMethod()
    hats[, "stage_2"]
  } else if(type == "stage1"){
    return(hats[, -ncol(hats)])
  } else {
    n <- model$nobs
    p <- model$p
    q <- model$q
    mean1 <- q/n
    mean2 <- p/n
    # hat2 <- lm.influence(model)$hat/mean2
    hat2 <- hats[, "stage_2"]/mean2
    # model[c("qr", "rank", "residuals", "coefficients")] <-
    #   list(model$qr1, model$rank1, na.omit(residuals(model, type="projected")),
    #        model$coefficients1)
    # hat1 <- lm.influence(model)$hat/mean1
    hat1 <- if (model$method == "OLS"){
      hats[, "stage_1"]
    } else {
      nendog <- ncol(hats) - 1
      if (type == "both"){
        exp(rowMeans(log(hats[, 1:nendog, drop=FALSE])))/mean1
      } else {
        apply(hats[, 1:nendog, drop=FALSE], 1, pmax)/mean1
      }
    }
    hat <- if (type == "both") {
      sqrt(hat1*hat2)
    } else {
      pmax(hat1, hat2)
    }
    mean2*hat
  }
  na.action <- model$na.action
  names(hatvalues) <- nms
  if (inherits(na.action, "exclude")) hatvalues[na.action]  <- NA
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

#' @rdname ivregDiagnostics
#' @method qqPlot influence.ivreg
#' @param distribution \code{"t"} (the default) or \code{"norm"}.
#' @param ylab The vertical axis label.
#' @export
qqPlot.influence.ivreg <- function(x,
                                   ylab = paste("Studentized Residuals(", deparse(substitute(x)), ")", sep = ""),
                                   distribution = c("t", "norm"), ...){
  distribution <- match.arg(distribution, c("t", "norm"))
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
}

#' @rdname ivregDiagnostics
#' @method model.matrix influence.ivreg
#' @export
model.matrix.influence.ivreg <- function(object, ...){
  object$model
}

#' @rdname ivregDiagnostics
#' @importFrom car avPlots
#' @export
avPlots.ivreg <- function(model, terms, ...){
  if (missing(terms)){
    terms <- formula(model, "regressors")[c(1, 3)]
    avPlots(model, terms=terms, ...)
    return(invisible(NULL))
  }
  .Class <- "lm"
  NextMethod()
}

#' @rdname ivregDiagnostics
#' @importFrom car avPlot
#' @export
avPlot.ivreg <- function(model, ...){
  xz <- model.matrix(model, component = "projected")
  if(is.null(model$x)) model$x <- list()
  model$x$regressors <- xz
  .Class <- "lm"
  NextMethod()
}

#' @rdname ivregDiagnostics
#' @importFrom car mcPlots
#' @export
mcPlots.ivreg <- function(model, terms, ...){
  if (missing(terms)){
    terms <- formula(model, "regressors")[c(1, 3)]
    mcPlots(model, terms=terms, ...)
    return(invisible(NULL))
  }
  .Class <- "lm"
  NextMethod()
}

#' @rdname ivregDiagnostics
#' @importFrom car mcPlot
#' @export
mcPlot.ivreg <- function(model, ...){
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
  method <- match.arg(method, "case")
  NextMethod()
}

#' @rdname ivregDiagnostics
#' @importFrom car crPlots
#' @export
crPlots.ivreg <- function(model, terms, ...){
  if (missing(terms)){
    terms <- formula(model, "regressors")[c(1, 3)]
    crPlots(model, terms=terms, ...)
    return(invisible(NULL))
  }
  .Class <- "lm"
  NextMethod()
}

#' @rdname ivregDiagnostics
#' @importFrom car crPlot
#' @export
crPlot.ivreg <- function(model, ...){
  model$contrasts <- model$contrasts$regressors
  .Class <- "lm"
  NextMethod()
}

#' @rdname ivregDiagnostics
#' @importFrom car ceresPlots
#' @export
ceresPlots.ivreg <- function(model, terms, ...){
  if (missing(terms)){
    terms <- formula(model, "regressors")[c(1, 3)]
    ceresPlots(model, terms=terms, ...)
    return(invisible(NULL))
  }
  .Class <- "lm"
  NextMethod()
}

#' @rdname ivregDiagnostics
#' @importFrom car ceresPlot
#' @export
ceresPlot.ivreg <- function(model, ...){
  model$contrasts <- model$contrasts$regressors
  model$formula <- formula(model, "regressors")
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
}

#' @rdname ivregDiagnostics
#' @importFrom car qqPlot
#' @export
qqPlot.ivreg <- function(x, distribution=c("t", "norm"), ...){
  distribution <- match.arg(distribution, c("t", "norm"))
  rstudent <- rstudent(x)
  df <- df.residual(x)
  qqPlot(rstudent, distribution=distribution, df=df, ylab="Studentized Residuals", ...)
}

#' @rdname ivregDiagnostics
#' @importFrom car outlierTest
#' @export
outlierTest.ivreg <- function(x, ...){
  .Class <- "lm"
  NextMethod()
}

#' @rdname ivregDiagnostics
#' @importFrom car influencePlot
#' @export
influencePlot.ivreg <- function(x, ...){
  .Class <- "lm"
  NextMethod()
}

#' @rdname ivregDiagnostics
#' @importFrom car spreadLevelPlot
#' @export
spreadLevelPlot.ivreg <- function(x, main="Spread-Level Plot", ...){
  .Class <- "lm"
  NextMethod()
}

#' @rdname ivregDiagnostics
#' @importFrom car ncvTest
#' @export
ncvTest.ivreg <- function(model, ...){
    .Class <- "lm"
    NextMethod()
}

#' @rdname ivregDiagnostics
#' @importFrom stats deviance
#' @export
deviance.ivreg <- function(object, ...){
  .Class <- "lm"
  NextMethod()
}

# #' @rdname ivregDiagnostics
# #' @export
# hatvalues.rivreg <- function(model, ...){
#   wts <- weights(model)
#   rwts <- model$rweights
#   nms <- rownames(rwts)
#   rwts <- rwts[, ncol(rwts)]
#   n <- length(rwts)
#   wts <- if (is.null(wts)) rwts else wts*rwts
#   model$weights <- wts
#   class(model) <- "ivreg"
#   hat <- hatvalues(model, ...)
#   result <- rep(0, n)
#   result[wts > 0] <- hat
#   names(result) <- nms
#   result
# }

#' @rdname ivregDiagnostics
#' @export
influence.rivreg <- function(model, ...){
  mwts <- wts <- weights(model)
  rwts <- model$rweights
  nms <- rownames(rwts)
  rwts <- rwts[, ncol(rwts)]
  n <- length(rwts)
  wts <- if (is.null(wts)) rwts else wts*rwts
  model$weights <- wts
  class(model) <- "ivreg"
  infl <- influence(model)
  if (is.null(mwts)) mwts <- 1
  infl$rstudent <- mwts*residuals(model)/(infl$sigma * sqrt(1 - infl$hatvalues))
  cookd <- infl$cookd
  cookd[is.nan(cookd)] <- 0
  infl$cookd <- cookd
  dffits <- infl$dffits
  dffits[is.nan(dffits)] <- 0
  infl$dffits <- dffits
  infl
}

hatvalues.rlm <- function(model, ...){
  # this unexported method is necessary because
  # the inherited lm method doesn't take
  # account of robustness weights
  wts <- weights(model)
  if (is.null(wts)) wts <- 1
  rwts <- model$w
  nms <- names(y <- residuals(model))
  n <- length(rwts)
  wts <-  wts*rwts
  if (any(wts == 0)){
    X <- model.matrix(model)
    # using residuals for y as a place holder
    # to get WLS hatvalues
    model <- lm(y ~ X - 1, weights = wts)
  }
  hat <- lm.influence(model, do.coef = FALSE)$hat
  # doesn't include hatvalues for cases with 0 weights
  result <- rep(0, n)
  result[wts > 0] <- hat
  names(result) <- nms
  result
}
