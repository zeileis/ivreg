na.remove <- function(x){
  # remove NAs from a vector preserving names and returning a vector
  if (!is.vector(x)) "x is not a vector"
  x <- na.omit(x)
  nms <- names(x)
  x <- as.vector(x)
  if (!is.null(nms)) names(x) <- nms
  x
}

#' 2SLS Regression
#'
#' @param y The response variable, a vector.
#' @param X The explanatory-variables model matrix.
#' @param Z The instrumental-variables model matrix.
#' @param wt An optional weight vector.
#' @param singular.ok Is perfect collinearity tolerable in either stage regression?
#' @param qr Include the QR decomposition in the returned object.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{\code{n}}{number of cases}
#'   \item{\code{p}}{number of columns of \code{X}}
#'   \item{\code{q}}{number of columns of \code{Z}}
#'   \item{\code{qr}}{QR decomposition for stage-2 regression}
#'   \item{\code{rank}}{rank of model matrix for stage-2 regression}
#'   \item{\code{qr.1}}{QR decomposition for stage-1 regression}
#'   \item{\code{rank.1}}{rank of \code{Z}}
#'   \item{\code{coefficients}}{estimated regression coefficients}
#'   \item{\code{coefficients.1}}{coefficient matrix from stage-1 regression}
#'   \item{\code{vcov}}{covariance matrix of estimated coefficients}
#'   \item{\code{df.residual}}{residual degrees of freedom for stage-2 regression}
#'   \item{\code{sigma}}{residual standard deviation}
#'   \item{\code{residuals}}{response residuals}
#'   \item{\code{fitted}}{model fitted values}
#'   \item{\code{residuals.1}}{matrix of stage-1 residuals}
#'   \item{\code{fitted.1}}{matrix of stage-1 fitted values}
#'   \item{\code{residuals.2}}{stage-2 residuals}
#'   \item{\code{fitted.2}}{stage-2 fitted values}
#'   \item{\code{weights}}{model weights or \code{NULL}}
#' }
#'
#' @importFrom stats lsfit naresid napredict
#' @export
#'
#' @description The work-horse 2SLS function, not normally meant to be called directly
#' but rather via \code{\link{ivreg2}}.
#'
#' @author John Fox \email{jfox@mcmaster.ca}
#'
#' @seealso \code{\link{ivreg2}}
#'
#' @examples
#' m <- with(Kmenta, fitivreg2(Q, cbind(1, P, D), cbind(1, D, F, A)))
#' names(m)
#' m$coefficients # estimates
#' sqrt(diag(m$vcov)) # standard errors
fitivreg2 <- function(y, X, Z, wt=NULL, singular.ok=FALSE, qr=TRUE){

  # handle NAs

  na.action <- which(rowSums(is.na(cbind(y, X, Z, wt))) > 0)
  if (length(na.action) > 0) {
    class(na.action) <- "exclude" # to be used for fits and residuals
    warning(length(na.action), " cases with missing values deleted")
    y <- y[-na.action]
    X <- X[-na.action, , drop=FALSE]
    Z <- Z[-na.action, , drop=FALSE]
    if (!is.null(wt)) wt <- wt[-na.action]
  }
  else na.action <- NULL

  rnames <- rownames(X)
  names(y) <- rnames

  # stage 1 regression

  stage1 <- suppressWarnings(lsfit(Z, X, wt=wt, intercept=FALSE))
  rk.1 <- stage1$qr$rank
  if (rk.1 < ncol(Z)){
    msg <- paste0("rank of stage-1 model matrix = ", rk.1,
                  " < number of stage-1 coefficients = ", ncol(Z))
    if (singular.ok) warning(msg) else stop(msg)
  }
  residuals.1 <- stage1$residuals
  rownames(residuals.1) <- rnames

  # stage 2 regression

  X.fit <- X - stage1$residuals
  colnames(X.fit) <- colnames(X)
  stage2 <- suppressWarnings(lsfit(X.fit, y, wt=wt, intercept=FALSE))
  rk.2 <- stage2$qr$rank
  if (rk.2 < ncol(X.fit)){
    msg <- paste0("rank of stage-2 model matrix = ", rk.2,
                  " < number of stage-2 coefficients = ", ncol(X.fit))
    if (singular.ok) warning(msg) else stop(msg)
  }
  residuals.2 <- stage2$residuals
  names(residuals.2) <- rnames

  fitted <- as.vector(X %*% stage2$coef)
  names(fitted) <- rnames
  residuals <- y - fitted
  names(residuals) <- rnames
  p <- ncol(X)
  n <- nrow(X)
  df.res <- n - p
  sigma2 <- (if (is.null(wt)) sum(residuals^2) else sum(wt*residuals^2))/df.res
  vcov <- sigma2*chol2inv(stage2$qr$qr)
  rownames(vcov) <- colnames(vcov) <- names(stage2$coef)

  list(
    n              = n,
    p              = p,
    q              = ncol(Z),
    qr             = if (qr) stage2$qr else NULL,
    rank           = rk.2,
    qr.1           = if (qr) stage1$qr else NULL,
    rank.1         = rk.1,
    coefficients   = stage2$coef,
    coefficients.1 = stage1$coef,
    vcov           = vcov,
    df.residual    = df.res,
    sigma          = sqrt(sigma2),
    residuals      = naresid(na.action, residuals),
    fitted         = napredict(na.action, fitted),
    residuals.1    = naresid(na.action, residuals.1),
    fitted.1       = napredict(na.action, X.fit),
    residuals.2    = naresid(na.action, residuals.2),
    fitted.2       = napredict(na.action, y - stage2$residuals),
    weights        = if (!is.null(wt)) naresid(na.action, wt) else NULL
  )
}

#' 2SLS Least Squares Estimation of a Linear Model
#'
#' @param formula A model \link[stats]{formula}, as for \code{\link[stats]{lm}},
#' specifying the response variable and the right-hand side of the model.
#' @param instruments A one-sided formula, specifying the instrumental variables.
#' @param data An optional data frame, list, or environment containing the variables in
#' the model formula and instrumental-variables formula.
#' @param subset See \code{\link[stats]{lm}}.
#' @param weights An optional vector of inverse-variance weights; if absent the weights
#' are set to 1 for all cases.
#' @param na.action See \code{\link[stats]{lm}}.
#' @param contrasts See \code{\link[stats]{lm}}.
#' @param singular.ok If \code{FALSE} (the default) perfect collinearity in either stage of
#' the 2SLS regression produces an error.
#' @param model Include the model frame in the returned object.
#' @param x Include the model matrices for the two stages in the returned object.
#' @param y Include the response vector in the returned object.
#' @param qr Include the QR decompositions for the two stages in the returned object.
#' @param ... Not used.
#'
#' @return An object of class \code{"ivreg2"}, with the following elements in addition to
#' those returned by \code{\link{fitivreg2}}:
#' \describe{
#'   \item{\code{response.name}}{a character string with the name of the response
#'   variable or possibly the character representation of an expression evaluating
#'   to the response.}
#'   \item{\code{formula}}{the model formula.}
#'   \item{\code{intruments}}{the one-sided formula specifying the instrumental variables.}
#'   \item{\code{model.matrix}}{the model matrix corresponding to the right-hand side of
#'   the model formula.}
#'   \item{\code{y}}{the response-variable vector.}
#'   \item{\code{model.matrix.instruments}}{the model matrix for the first-stage regression.}
#'   \item{\code{wt.var}}{the name of the weight variable (for a weighted fit).}
#'   \item{\code{na.action}}{the na.action argument.}
#'   \item{\code{call}}{the function call.}
#'   \item{\code{contasts}}{contrasts for the model matrix.}
#'   \item{\code{contrasts.instruments}}{contrasts for the instrumental-variables
#'   model matrix.}
#'   \item{\code{xlevels}}{factor levels (if any) for the model matrix.}
#'   \item{\code{xlevels.instruments}}{factor levels (if any) for the
#'   instrumental variables model matrix.}
#'   \item{\code{terms}}{the terms object for the model matrix.}
#'   \item{\code{terms.instruments}}{the terms object for the instrumental
#'   variables model matrix.}
#'   \item{\code{model}}{the model frame.}
#'   }
#'
#' @description   Provides a formula-based interface for 2SLS estimation of
#' a linear model. Computations are done by \code{\link{fitivreg2}}. The returned object has
#' the necessary information for computing a variety of
#' \link[=2SLS_Diagnostics]{regression diagnostics}.
#'
#' @author John Fox \email{jfox@mcmaster.ca}
#'
#' @seealso \code{\link{fitivreg2}}, \code{\link[stats]{lm}}, \code{\link[stats]{formula}}, \code{\link{2SLS_Methods}},
#' \code{\link{2SLS_Diagnostics}}
#'
#' @importFrom stats model.weights .getXlevels
#' @export
#'
#' @examples
#'summary(ivreg2(Q ~ P + D, ~ D + F + A, data=Kmenta))     # demand equation
#'summary(ivreg2(Q ~ P + F + A, ~ D + F + A, data=Kmenta)) # supply equation
#'
ivreg2 <- function (formula, instruments=rhs(formula), data, subset, weights,
                    na.action=getOption("na.action"), contrasts = NULL,
                    singular.ok=FALSE, model=TRUE, x=TRUE, y=TRUE, qr=TRUE, ...){
  rhs <- function(formula) if (length(formula) == 3) formula[-2] else formula
  combineFormulas <- function(formula1, formula2){
    rhs <- as.character(formula2)[length(formula2)]
    formula2 <- paste("~ . +", rhs)
    update(formula1, formula2)
  }
  cl <- match.call()
  m <- match.call(expand.dots = FALSE)
  if (is.matrix(eval(m$data, parent.frame()))) m$data <- as.data.frame(data)
  response.name <- deparse(formula[[2]])
  m$formula <- combineFormulas(formula, instruments)
  m$instruments <- m$contrasts <- m$singular.ok <- NULL
  m[[1]] <- as.name("model.frame")
  mt.model <- mt.instruments <- m
  mt.model$formula <- formula
  mt.instruments$formula <- instruments
  mf <- eval(m, parent.frame())
  mt.model <- attr(eval(mt.model, parent.frame()), "terms")
  mt.instruments <- attr(eval(mt.instruments, parent.frame()), "terms")
  na.act <- attr(mf, "na.action")
  w <- as.vector(model.weights(mf))
  wt.var <- if (!is.null(w)) deparse(substitute(weights)) else NULL
  Z <- model.matrix(instruments, data = mf, contrasts)
  y. <- mf[, response.name]
  X <- model.matrix(formula, data = mf, contrasts)
  names(y.) <- rownames(X)
  result <- fitivreg2(y., X, Z, w, singular.ok=singular.ok, qr=qr)
  result <- c(result, list(
    response.name            = response.name,
    formula                  = formula,
    instruments              = instruments,
    model.matrix             = if (x) X else NULL,
    y                        = if (y) y. else NULL,
    model.matrix.instruments = if (x) Z else NULL,
    wt.var                   = wt.var,
    na.action                = na.act,
    call                     = cl,
    contrasts                = attr(X, "contrasts"),
    contrasts.instruments    = attr(Z, "contrasts"),
    xlevels                  = .getXlevels(mt.model, mf),
    xlevels.instruments      = .getXlevels(mt.instruments, mf),
    terms                    = mt.model,
    terms.instruments        = mt.instruments,
    model                    = if (model) mf else NULL
  ))
  class(result) <- c("ivreg2", "lm")
  result
}


#' Methods for \code{"ivreg2"} Objects
#' @aliases 2SLS_Methods
#' @description Various methods for processing \code{"ivreg2"} objects; for diagnostic methods,
#'   see \link{2SLS_Diagnostics}.
#' @param object An object of class \code{"ivreg2"}.
#' @param type Type of object desired, varies by method:
#' \describe{
#'   \item{\code{model.matrix}}{\code{"model"} (the default), \code{"instruments"}, or
#'     \code{"stage2"}.}
#'   \item{\code{residuals}}{\code{"response"}, \code{"stage1"}, \code{"stage2"},
#'     \code{"working"} (equivalent to \code{"response"}), \code{"deviance"},
#'     \code{"pearson"} (equivalent to \code{"deviance"}), or \code{"partial"}.}
#'   \item{\code{fitted}}{\code{"model"}, \code{"stage1"}, or \code{"stage2"}.}
#'   }
#' @param ... to match generics, not generally used.
#'
#' @importFrom stats model.matrix
#' @export
#' @method model.matrix ivreg2
#' @seealso \code{\link{ivreg2}}, \link{2SLS_Diagnostics},
#'   \code{\link[sandwich]{sandwich}}
#' @examples
#' kmenta.eq1 <- ivreg2(Q ~ P + D, ~ D + F + A, data=Kmenta)
#' coef(kmenta.eq1) # estimates
#' sqrt(diag(vcov(kmenta.eq1))) # std. errors
#' summary(kmenta.eq1)
#' summary(kmenta.eq1, vcov.=sandwich::sandwich) # sandwich SEs
#' plot(fitted(kmenta.eq1), residuals(kmenta.eq1)) # residuals vs fitted values
model.matrix.ivreg2 <- function(object, type=c("model", "instruments", "stage2"), ...){
  type <- match.arg(type)
  switch(type,
         model = object$model.matrix,
         instruments = object$model.matrix.instruments,
         stage2 = object$fitted.1)
}

#' @rdname model.matrix.ivreg2
#' @param model An object of class \code{"ivreg2"} or \code{"influence.ivreg2"}.
#' @importFrom car avPlot
#' @export
avPlot.ivreg2 <- function(model, ...){
  model$model.matrix <- model.matrix(model, type="stage2")
  NextMethod()
}

#' @rdname model.matrix.ivreg2
#' @importFrom stats vcov
#' @export
vcov.ivreg2 <- function(object, ...) {
  object$vcov
}

#' @rdname model.matrix.ivreg2
#' @importFrom stats residuals predict
#' @export
residuals.ivreg2 <- function(object, type=c("response", "stage1", "stage2", "working",
                                          "deviance", "pearson", "partial"), ...){
  type <- match.arg(type)
  w <- object$weights
  if (is.null(w)) w <- 1
  res <- switch(type,
         working  =,
         response  = object$residuals,
         deviance =,
         pearson  = sqrt(w)*object$residuals,
         stage1   = object$residuals.1,
         stage2   = object$residuals.2,
         partial  = object$residuals + predict(object, type = "terms"))
  naresid(object$na.action, res)
}

#' @rdname model.matrix.ivreg2
#' @importFrom stats fitted
#' @export
fitted.ivreg2 <- function(object, type=c("model", "stage1", "stage2"), ...){
  type <- match.arg(type)
  switch(type,
         model  = napredict(object$na.action, object$fitted),
         stage1 = napredict(object$na.action, object$fitted.1),
         stage2 = napredict(object$na.action, object$fitted.2))
}

#' @rdname model.matrix.ivreg2
#' @param x An object of class \code{"ivreg2"}.
#' @param digits Digits to print.
#' @importFrom stats coef
#' @export
print.ivreg2 <- function (x, digits = getOption("digits") - 2, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\n")
  print(coef(x), digits=digits)
  cat("\nResidual standard deviation =",
      paste0(format(x$sigma, digits=digits), ",  R-squared analog ="), Rsq(x))
  invisible(x)
}

#' @rdname model.matrix.ivreg2
#' @param vcov. Function to compute the coefficient covariance matrix or the matrix itself;
#'   the default is the function \code{vcov}; set, e.g., to \code{\link[sandwich]{sandwich}}
#'   (from the \pkg{sandwich} package) to get robust coefficient standard errors.
#'
#' @importFrom stats pt getCall formula
#' @export
summary.ivreg2 <- function (object, digits = getOption("digits") - 2, vcov.=vcov, ...){
  df <- df.residual(object)
  std.errors <- if (is.function(vcov.)) {
    sqrt(diag(V <- vcov.(object)))
  } else {
    if (!is.matrix(vcov.)) stop("vcov. must be a function or a matrix")
    if (!(ncol(vcov.) == nrow(vcov.)) || ncol(vcov.) != length(coef(object)))
      stop("vcov. matrix is of the wrong dimensions")
    sqrt(diag(V <- vcov.))
  }
  b <- coef(object)
  t <- b/std.errors
  p <- 2 * (1 - pt(abs(t), df))
  table <- cbind(b, std.errors, t, p)
  rownames(table) <- names(b)
  colnames(table) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  intercept <- which(names(b) == "(Intercept)")
  b <- b[-intercept]
  V <- V[-intercept, -intercept]
  F <-  (b %*% solve(V) %*% b) / length(b)
  pval <- pf(F, length(b), df, lower.tail=FALSE)
  result <- list(call=object$call,
                 residuals = summary(residuals(object)),
                 coefficients = table, digits = digits, sigma = object$s,
                 df = df, dfn=length(b), na.action=object$na.action, r2=Rsq(object),
                 r2adj=Rsq(object, adjust=TRUE), F=F, pval=pval)
  class(result) <- "summary.ivreg2"
  result
}

#' @rdname model.matrix.ivreg2
#' @method print summary.ivreg2
#' @importFrom stats printCoefmat
#' @export
print.summary.ivreg2 <- function (x, ...) {
  digits <- x$digits
  cat("Call:\n")
  print(x$call)
  cat("\nResiduals:\n")
  print(round(x$residuals, digits))
  cat("\nCoefficients:\n")
  printCoefmat(x$coefficients, digits = digits)
  cat("\nResidual standard deviation =", round(x$sigma, digits),
            "on", x$df, "degrees of freedom")
  cat("\nR-squared analog =", paste0(format(x$r2, digits=digits),
                                     ",  Adjusted R-squared ="),
      format(x$r2adj, digits=digits))
  cat("\nF =", paste0(format(x$F, digits=digits), " on ", x$dfn, " and ", x$df),
      "degrees of freedom, p", if (x$pval < .Machine$double.eps) "" else "=",
      format.pval(x$pval, digits=digits), "\n\n")
  if (!is.null(x$na.action)){
    cat(paste(length(x$na.action), "cases deleted due to NAs\n\n"))
  }
  invisible(x)
}

#' @rdname model.matrix.ivreg2
#' @param model.2 A second object of class \code{"ivreg2"}.
#' @param s2 the estimated error variance (optional);
#'   if not specified, taken from the larger model.
#' @param dfe the degrees of freedom for error (optional);
#'   if not specified, taken from the larger model.
#' @importFrom stats anova pf
#' @export
anova.ivreg2 <- function(object, model.2, s2, dfe, ...){
  if (!inherits(model.2, "ivreg2")) stop('requires two models of class ivreg2')
  s2.1 <- object$s^2
  dfe.1 <- df.residual(object)
  s2.2 <- model.2$s^2
  dfe.2 <- df.residual(model.2)
  SS.1 <- s2.1 * dfe.1
  SS.2 <- s2.2 * dfe.2
  SS <- abs(SS.1 - SS.2)
  Df <- abs(dfe.2 - dfe.1)
  if (missing(s2)){
    s2 <- if (dfe.1 > dfe.2) s2.2 else s2.1
    f <- (SS/Df) / s2
    RSS <- c(SS.1, SS.2)
    Res.Df <- c(dfe.1, dfe.2)
    SS <- c(NA, SS)
    P <- c(NA, 1 - pf(f, Df, min(dfe.1, dfe.2)))
    Df <- c(NA, Df)
    f <- c(NA, f)
    rows <- c("Model 1", "Model 2")
  }
  else{
    f <- (SS/Df) / s2
    RSS <- c(SS.1, SS.2, s2*dfe)
    Res.Df <- c(dfe.1, dfe.2, dfe)
    SS <- c(NA, SS, NA)
    P <- c(NA, 1 - pf(f, Df, dfe), NA)
    Df <- c(NA, Df, NA)
    f <- c(NA, f, NA)
    rows <- c("Model 1", "Model 2", "Error")
  }
  table <- data.frame(Res.Df, RSS, Df, SS, f, P)
  head.1 <- paste("Model 1: ",format(object$formula), "  Instruments:",
                  format(object$instruments))
  head.2 <- paste("Model 2: ",format(model.2$formula), "  Instruments:",
                  format(model.2$instruments))
  names(table) <- c("Res.Df", "RSS", "Df", "Sum of Sq", "F", "Pr(>F)")
  row.names(table) <- rows
  structure(table, heading = c("Analysis of Variance", "", head.1, head.2, ""),
            class = c("anova", "data.frame"))
}

#' @rdname model.matrix.ivreg2
#' @param formula. Updated model formula.
#' @param instruments. Updated one-sided formula for the instrumental variables.
#' @param evaluate If \code{TRUE} (the default) evaluate the updated model;
#'   if \code{FALSE} simply generate the updated call.
#' @importFrom stats update
#' @export
update.ivreg2 <- function (object, formula., instruments., ..., evaluate=TRUE){
  # adapted from stats::update.default()
  if (is.null(call <- getCall(object)))
    stop("need an object with call component")
  extras <- match.call(expand.dots = FALSE)$...
  if (!missing(formula.))
    call$formula <- update(formula(object), formula.)
  if (!missing(instruments.))
    call$instruments <- update(object$instruments, instruments.)
  if (length(extras)) {
    existing <- !is.na(match(names(extras), names(call)))
    for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
    if (any(!existing)) {
      call <- c(as.list(call), extras[!existing])
      call <- as.call(call)
    }
  }
  if (evaluate)
    eval(call, parent.frame())
  else call
}

#' @rdname model.matrix.ivreg2
#' @importFrom sandwich bread
#' @export
bread.ivreg2 <- function(x, ...){
  x$vcov*x$n/x$sigma^2
}

diagprod <- function(d, X){
  # equivalent to diag(d) %*% X
  if (!is.vector(d)) stop("d is not a vector")
  if (!is.matrix(X)) stop("X is not a matrix")
  if (length(d) != nrow(X)) stop("d and X not conformable")
  d*X
}

#' @rdname model.matrix.ivreg2
#' @importFrom sandwich estfun
#' @export
estfun.ivreg2 <- function (x, ...) {
  if (x$rank < x$p) stop("second stage model matrix is of deficient rank")
  w <- weights(x)
  if (is.null(w)) w <- 1
  diagprod(na.remove(w*residuals(x)), model.matrix(x, type="stage2"))
}

#' R-Squares Analog
#'
#' @param model A model object.
#' @param ... Possible arguments for specific methods.
#' @importFrom stats na.omit model.response model.frame df.residual weights
#' @export
#'
#' @examples
#' Rsq(ivreg2(Q ~ P + D, ~ D + F + A, data=Kmenta))
Rsq <- function(model, ...){
  UseMethod("Rsq")
}

#' @rdname Rsq
#' @export
#' @param adjusted If \code{TRUE} (the default is \code{FALSE}) return the \eqn{R^2} adjusted
#' for degrees of freedom.
Rsq.default <- function(model, adjusted=FALSE, ...){
  w <- weights(model)
  if (is.null(w)) w <- 1
  SSE <- sum(w*residuals(model)^2, na.rm=TRUE)
  y <- model$y
  SST <- sum(na.omit(w)*(y - mean(y))^2)
  if (adjusted) {
    1 - (SSE/df.residual(model))/(SST/(model$n - 1))
  }
  else {
    1 - SSE/SST
  }
}
