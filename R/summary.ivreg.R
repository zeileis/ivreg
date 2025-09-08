#' Summary and Inference Methods for \code{"ivreg"} Objects
#' @aliases summary.ivreg print.summary.ivreg anova.ivreg confint.ivreg Anova.ivreg linearHypothesis.ivreg
#' @description Summary method, including Wald tests and (by default) certain diagnostic tests, for
#' \code{"ivreg"} model objects, as well as other related inference functions.
#' @seealso \code{\link{ivreg}}, \code{\link{ivreg.fit}}, \code{\link{ivregDiagnostics}}
#' @param object,object2,model,mod An object of class \code{"ivreg"}.
#' @param x An object of class \code{"summary.ivreg"}.
#' @param component Character indicating \code{"stage2"} or \code{"stage1"}.
#' @param digits Minimal number of significant digits for printing.
#' @param signif.stars Show "significance stars" in summary output?
#' @param vcov. Optionally either a coefficient covariance matrix or a function to compute such a covariance
#'   matrix from fitted \code{ivreg} model objects. If \code{NULL} (the default) the standard covariance matrix
#'   (based on the information matrix) is used. Alternatively, covariance matrices (e.g., clustered and/or
#'   heteroscedasticity-consistent) can be plugged in to adjust Wald tests or confidence intervals etc.
#'   In \code{summary}, if \code{diagnostics = TRUE}, \code{vcov.} must be a function (not a matrix) because
#'   the alternative covariances are also needed for certain auxiliary models in the diagnostic tests.
#'   If \code{vcov.} is a function, the \code{...} argument can be used to pass on further arguments to
#'   this function.
#' @param df For \code{summary}, optional residual degrees of freedom to use in computing model summary. 
#' @param diagnostics Report 2SLS "diagnostic" tests in model summary (default is \code{TRUE}). 
#' These tests are not to be confused with the \emph{regression diagnostics} provided elsewhere in the \pkg{ivreg}
#' package: see \code{\link{ivregDiagnostics}}.
#' @param test,test.statistic Test statistics for ANOVA table computed by \code{anova}, \code{\link[car]{Anova}},
#' or \code{\link[car]{linearHypothesis}}. Only \code{test = "F"} is supported by \code{anova}; this is also
#' the default for \code{Anova} and \code{linearHypothesis}, which also allow \code{test = "Chisq"} for
#' asymptotic tests.
#' @param hypothesis.matrix,rhs For formulating a linear hypothesis; see the documentation
#' for \code{\link[car]{linearHypothesis}} for details.
#' @param complete If \code{TRUE}, the default, the returned coefficient vector (for \code{coef}) or coefficient-covariance matrix (for \code{vcov}) includes elements for aliased regressors.
#' @param parm  parameters for which confidence intervals are to be computed; a vector or numbers or names; the default is all parameters.
#' @param level confidence level; the default is \code{0.95}.
#' @param ... arguments to pass down.
#' @examples
#' \dontshow{ if(!requireNamespace("sandwich")) {
#'   if(interactive() || is.na(Sys.getenv("_R_CHECK_PACKAGE_NAME_", NA))) {
#'     stop("not all packages required for the example are installed")
#'   } else q() }}
#' ## data and model
#' data("CigaretteDemand", package = "ivreg")
#' m <- ivreg(log(packs) ~ log(rincome) | log(rprice) | salestax, data = CigaretteDemand)
#' 
#' ## summary including diagnostics
#' summary(m)
#' 
#' ## replicate global F test from summary (against null model) "by hand"
#' m0 <- ivreg(log(packs) ~ 1, data = CigaretteDemand)
#' anova(m0, m)
#' 
#' ## or via linear hypothesis test
#' car::linearHypothesis(m, c("log(rincome)", "log(rprice)"))
#'
#' ## confidence intervals
#' confint(m)
#'
#' ## just the Wald tests for the coefficients
#' library("lmtest")
#' coeftest(m)
#' 
#' ## plug in a heteroscedasticity-consistent HC1 covariance matrix (from sandwich)
#' library("sandwich")
#' ## - as a function passing additional type argument through ...
#' coeftest(m, vcov = vcovHC, type = "HC1")
#' ## - as a function without additional arguments
#' hc1 <- function(object, ...) vcovHC(object, type = "HC1", ...)
#' coeftest(m, vcov = hc1)
#' ## - as a matrix
#' vc1 <- vcovHC(m, type = "HC1")
#' coeftest(m, vcov = vc1)
#' 
#' ## in summary() with diagnostics = TRUE use one of the function specifications,
#' ## the matrix is only possible when diagnostics = FALSE
#' summary(m, vcov = vcovHC, type = "HC1")     ## function + ...
#' summary(m, vcov = hc1)                      ## function
#' summary(m, vcov = vc1, diagnostics = FALSE) ## matrix
#'
#' ## in confint() and anova() any of the three specifications can be used
#' anova(m0, m, vcov = vcovHC, type = "HC1")   ## function + ...
#' anova(m0, m, vcov = hc1)                    ## function
#' anova(m0, m, vcov = vc1)                    ## matrix
#'
#' @importFrom stats model.matrix vcov .vcov.aliased anova quantile weighted.mean lm lm.fit lm.wfit pchisq 
#' @importFrom lmtest coeftest waldtest waldtest.default lrtest lrtest.default
#' @importFrom car Anova linearHypothesis

#' @rdname summary.ivreg
#' @importFrom stats qnorm qt
#' @export
confint.ivreg <- function (object, parm, level = 0.95,
  component = c("stage2", "stage1"), complete = TRUE, vcov. = NULL, df = NULL, ...) 
{
  component <- match.arg(component, c("stage2", "stage1"))
  est <- coef(object, component = component, complete = complete)
  se <- if(is.null(vcov.)) {
    vcov(object, component = component, complete = complete)
  } else {
    if(is.function(vcov.)) {
      vcov.(object, ...)
    } else {
      vcov.
    }
  }
  se <- sqrt(diag(se))
  a <- (1 - level)/2
  a <- c(a, 1 - a)
  if(is.null(df)) df <- if(component == "stage2") object$df.residual else object$df.residual1
  crit <- if(is.finite(df) && df > 0) qt(a, df = df) else qnorm(a)
  ci <- cbind(est + crit[1L] * se, est + crit[2L] * se)
  colnames(ci) <- paste(format(100 * a, trim = TRUE, scientific = FALSE, digits = 3L), "%")
  if(missing(parm) || is.null(parm)) parm <- seq_along(est)
  if(is.character(parm)) parm <- which(names(est) %in% parm)
  ci <- ci[parm, , drop = FALSE]
  ci
}

#' @rdname summary.ivreg
#' @export
summary.ivreg <- function(object, vcov. = NULL, df = NULL, diagnostics = NULL, ...)
{
  # prevent some "inherited" "lm" methods from failing
  if(is.null(diagnostics)) diagnostics <- (length(object$endogenous) > 0L) && (length(object$instruments) > 0L)
  if(diagnostics && ((length(object$endogenous) <= 0L) || (length(object$instruments) <= 0L) || (length(formula(object, component="instruments")) <= 0L))) {
    diagnostics <- FALSE
    warning("diagnostics cannot be computed without endogenous/instrument variables")
  }
  
  ## weighted residuals
  res <- object$residuals
  y <- object$fitted.values + res
  n <- NROW(res)
  w <- object$weights
  if(is.null(w)) w <- rep(1, n)
  res <- res * sqrt(w)

  ## R-squared
  rss <- sum(res^2)
  if(attr(object$terms$regressors, "intercept")) {
    tss <- sum(w * (y - weighted.mean(y, w))^2)
    dfi <- 1    
  } else {
    tss <- sum(w * y^2)
    dfi <- 0
  }
  r.squared <- 1 - rss/tss
  adj.r.squared <- 1 - (1 - r.squared) * ((n - dfi)/object$df.residual)
  
  ## degrees of freedom (for z vs. t test)
  if(is.null(df)) df <- object$df.residual
  if(!is.finite(df)) df <- 0
  if(df > 0 & (df != object$df.residual)) {
    df <- object$df.residual
  }

  ## covariance matrix processing
  ## - sanity check
  if(!is.null(vcov.) && !is.function(vcov.) && diagnostics) {
    warning("'vcov.' must be a function for extracting covariance matrix estimates when using 'diagnostics = TRUE', using 'vcov. = NULL' instead")
    vcov. <- NULL
  }
  ## - covariance matrix for model object
  if(is.null(vcov.)) 
    vc <- vcov(object)
  else {
    vc <- if(is.function(vcov.)) vcov.(object, ...) else vcov.
  }
  ## - extractor function (if available)
  vcfun <- if(!is.null(vcov.) && is.function(vcov.)) function(object) vcov.(object, ...) else NULL
  
  ## Wald test of each coefficient
  cf <- lmtest::coeftest(object, vcov. = vc, df = df, ...)
  attr(cf, "method") <- NULL
  class(cf) <- "matrix"
  
  ## Wald test of all coefficients
  Rmat <- if(attr(object$terms$regressors, "intercept"))
    cbind(0, diag(length(na.omit(coef(object))) - 1)) else diag(length(na.omit(coef(object))))
  waldtest <- car::linearHypothesis(object, Rmat, vcov. = vc, test = ifelse(df > 0, "F", "Chisq"), singular.ok = TRUE)
  waldtest <- c(waldtest[2, "F"], waldtest[2, "Pr(>F)"], waldtest[2, "Df"], if(df > 0) waldtest[2, "Res.Df"] else NULL)
  
  ## diagnostic tests
  diag <- if(diagnostics) ivdiag(object, vcov. = vcfun) else NULL
  
  rval <- list(
    call = object$call,
    terms = object$terms,
    residuals = res,
    weights = object$weights,
    coefficients = cf,
    sigma = object$sigma,
    df = c(object$rank, if(df > 0) df else Inf, object$rank), ## aliasing
    r.squared = r.squared,
    adj.r.squared = adj.r.squared,
    waldtest = waldtest,
    vcov = vc,
    diagnostics = diag)
    
  class(rval) <- "summary.ivreg"
  return(rval)
}

#' @rdname summary.ivreg
#' @importFrom stats printCoefmat
#' @export
#' @method print summary.ivreg
print.summary.ivreg <- function(x, digits = max(3, getOption("digits") - 3), 
    signif.stars = getOption("show.signif.stars"), ...)
{
  cat("\nCall:\n")
  cat(paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  cat(if(!is.null(x$weights) && diff(range(x$weights))) "Weighted ", "Residuals:\n", sep = "")      
  if(NROW(x$residuals) > 5L) {
      nam <- c("Min", "1Q", "Median", "3Q", "Max")
      rq <- if(length(dim(x$residuals)) == 2) 
	  structure(apply(t(x$residuals), 1, quantile), dimnames = list(nam, dimnames(x$residuals)[[2]]))
      else structure(quantile(x$residuals), names = nam)
      print(rq, digits = digits, ...)
  } else {
      print(x$residuals, digits = digits, ...)
  }

  cat("\nCoefficients:\n")
  printCoefmat(x$coefficients, digits = digits, signif.stars = signif.stars,
    signif.legend = signif.stars & is.null(x$diagnostics), na.print = "NA", ...)

  if(!is.null(x$diagnostics)) {
    cat("\nDiagnostic tests:\n")
    printCoefmat(x$diagnostics, cs.ind = 1L:2L, tst.ind = 3L,
      has.Pvalue = TRUE, P.values = TRUE, digits = digits,
      signif.stars = signif.stars, na.print = "NA", ...)  
  }

  cat("\nResidual standard error:", format(signif(x$sigma, digits)),
    "on", x$df[2L], "degrees of freedom\n")

  cat("Multiple R-Squared:", formatC(x$r.squared, digits = digits))
  cat(",\tAdjusted R-squared:", formatC(x$adj.r.squared, digits = digits),
    "\nWald test:", formatC(x$waldtest[1L], digits = digits),
    "on", x$waldtest[3L], if(length(x$waldtest) > 3L) c("and", x$waldtest[4L]) else NULL,       
    "DF,  p-value:", format.pval(x$waldtest[2L], digits = digits), "\n\n")

  invisible(x)
}

#' @rdname summary.ivreg
#' @export    
anova.ivreg <- function(object, object2, test = "F", vcov. = NULL, ...)
{
  ## pass '...' (if any) to vcov. (if it is a function), ignore otherwise
  dots <- list(...)
  vc <- if(length(list(...)) > 0L && is.function(vcov.)) {
    function(object, ...) do.call("vcov.", c(list(object), dots, list(...)))
  } else {
    vcov.
  }
  
  ## call lmtest::waldtest
  rval <- waldtest(object, object2, test = test, vcov = vc)

  ## format as "anova" object if default vcov is used
  if(is.null(vc)) {
    head <- attr(rval, "heading")
    head[1] <- "Analysis of Variance Table\n"
    rss <- sapply(list(object, object2), function(x) sum(residuals(x)^2))
    dss <- c(NA, -diff(rss))
    rval <- cbind(rval, cbind("RSS" = rss, "Sum of Sq" = dss))[,c(1L, 5L, 2L, 6L, 3L:4L)]
    attr(rval, "heading") <- head
    class(rval) <- c("anova", "data.frame")
  }
  return(rval)
}

#' @rdname summary.ivreg
#' @importFrom car Anova
#' @export
Anova.ivreg <- function(mod, test.statistic = c("F", "Chisq"), vcov. = NULL, ...) {
  test.statistic <- match.arg(test.statistic, c("F", "Chisq"))
  NextMethod(test.statistic = test.statistic)
}

#' @rdname summary.ivreg
#' @export
linearHypothesis.ivreg <- function(model, hypothesis.matrix,
  rhs = NULL, test = c("F", "Chisq"), vcov. = NULL, ...)
{
  test <- match.arg(test, c("F", "Chisq"))
  NextMethod(test = test)
}

#' @importFrom stats model.frame model.response pf
ivdiag <- function(obj, vcov. = NULL) {
  ## extract data
  y <- model.response(model.frame(obj))
  x <- model.matrix(obj, component = "regressors")
  z <- model.matrix(obj, component = "instruments")
  w <- weights(obj)
  
  ## names of "regressors" and "instruments"
  xnam <- colnames(x)
  znam <- colnames(z)

  ## endogenous/instrument variables
  endo <- obj$endogenous
  inst <- obj$instruments
  if((length(endo) <= 0L) | (length(inst) <= 0L))
    stop("no endogenous/instrument variables")

  ## return value
  rval <- matrix(NA, nrow = length(endo) + 2L, ncol = 4L)
  colnames(rval) <- c("df1", "df2", "statistic", "p-value")
  rownames(rval) <- c(if(length(endo) > 1L) paste0("Weak instruments (", xnam[endo], ")") else "Weak instruments",
    "Wu-Hausman", "Sargan")
  
  ## convenience functions
  lmfit <- function(x, y, w = NULL) {
    rval <- if(is.null(w)) lm.fit(x, y) else lm.wfit(x, y, w)
    rval$x <- x
    rval$y <- y
    return(rval)
  }
  rss <- function(obj, weights = NULL) if(is.null(weights)) sum(obj$residuals^2) else sum(weights * obj$residuals^2)
  wald <- function(obj0, obj1, vcov. = NULL, weights = NULL) {
    df <- c(obj1$rank - obj0$rank, obj1$df.residual)
    if(!is.function(vcov.)) {
      w <- ((rss(obj0, w) - rss(obj1, w)) / df[1L]) / (rss(obj1, w)/df[2L])
    } else {
      if(NCOL(obj0$coefficients) > 1L) {
        cf0 <- structure(as.vector(obj0$coefficients),
	  .Names = c(outer(rownames(obj0$coefficients), colnames(obj0$coefficients), paste, sep = ":")))
        cf1 <- structure(as.vector(obj1$coefficients),
	  .Names = c(outer(rownames(obj1$coefficients), colnames(obj1$coefficients), paste, sep = ":")))
      } else {
        cf0 <- obj0$coefficients
        cf1 <- obj1$coefficients
      }
      cf0 <- na.omit(cf0)
      cf1 <- na.omit(cf1)
      ovar <- which(!(names(cf1) %in% names(cf0)))
      vc <- vcov.(lm(obj1$y ~ 0 + obj1$x, weights = w))
      w <- t(cf1[ovar]) %*% solve(vc[ovar,ovar]) %*% cf1[ovar]
      w <- w / df[1L]
    }
    pval <- pf(w, df[1L], df[2L], lower.tail = FALSE)
    c(df, w, pval)
  }
    
  # Test for weak instruments
  for(i in seq_along(endo)) {
    aux0 <- lmfit(z[, -inst, drop = FALSE], x[, endo[i]], w)
    aux1 <- lmfit(z,                        x[, endo[i]], w)
    rval[i, ] <- wald(aux0, aux1, vcov. = vcov., weights = w)
  }

  ## Wu-Hausman test for endogeneity
  if(length(endo) > 1L) aux1 <- lmfit(z, x[, endo], w)
  xfit <- as.matrix(aux1$fitted.values)
  colnames(xfit) <- paste("fit", colnames(xfit), sep = "_")
  auxo <- lmfit(      x,        y, w)
  auxe <- lmfit(cbind(x, xfit), y, w)
  rval[nrow(rval) - 1L, ] <- wald(auxo, auxe, vcov. = vcov., weights = w)

  ## Sargan test of overidentifying restrictions 
  r <- residuals(obj)  
  auxs <- lmfit(z, r, w)
  rssr <- if(is.null(w)) sum((r - mean(r))^2) else sum(w * (r - weighted.mean(r, w))^2)
  rval[nrow(rval), 1L] <- length(inst) - length(endo)
  if(rval[nrow(rval), 1L] > 0L) {
    rval[nrow(rval), 3L] <- length(r) * (1 - rss(auxs, w)/rssr)
    rval[nrow(rval), 4L] <- pchisq(rval[nrow(rval), 3L], rval[nrow(rval), 1L], lower.tail = FALSE)
  }

  return(rval)
}
