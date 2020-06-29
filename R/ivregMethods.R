#' Methods for \code{"ivreg"} Objects
#' @aliases ivreg_Methods vcov.ivreg bread.ivreg estfun.ivreg terms.ivreg model.matrix.ivreg predict.ivreg
#' print.ivreg summary.ivreg print.summary.ivreg anova.ivreg update.ivreg residuals.ivreg Effect.ivreg 
#' formula.ivreg
#' @description Various methods for processing \code{"ivreg"} objects; for diagnostic methods,
#'   see \code{\link{ivregDiagnostics}}.
#' @seealso \code{\link{ivreg}}, \code{\link{ivreg.fit}}, \code{\link{ivregDiagnostics}}
#' @param object,object2,model,mod An object of class \code{"ivreg"}.
#' @param x An object of class \code{"ivreg"} or \code{"summary.ivreg"}.
#' @param component For \code{\link{terms}}, \code{"regressors"}, \code{"instruments"}, or \code{"full"}; 
#' for \code{\link{model.matrix}}, \code{"projected"}, \code{"regressors"}, or \code{"instruments"};
#' for \code{\link{formula}}, \code{"regressors"}, \code{"instruments"},  or \code{"complete"}.
#' @param newdata Values of predictors for which to obtain predicted values.
#' @param na.action \code{na} method to apply to predictor values for predictions; default is \code{\link{na.pass}}.
#' @param digits For printing.
#' @param signif.stars Show "significance stars" in summary output.
#' @param vcov. Optional coefficient covariance matrix, or a function to compute the covariance matrix, to use in computing the model summary.
#' @param df Optional residual degrees of freedom to use in computing model summary.
#' @param diagnostics Report 2SLS "diagnostic" tests in model summary (default is \code{TRUE}). 
#' These tests are not to be confused with the \emph{regression diagnostics} provided elsewhere in the \pkg{ivreg}
#' package: see \code{\link{ivregDiagnostics}}.
#' @param test,test.statistic Test statistics for ANOVA table computed by \code{anova()}, \code{Anova()},
#' or \code{linearHypothesis()}. Only \code{test = "F"} is supported by \code{anova()}; this is also
#' the default for \code{Anova()} and \code{linearHypothesis()}, which also allow \code{test = "Chisq"} for
#' asymptotic tests.
#' @param formula. To update model.
#' @param evaluate If \code{TRUE}, the default, the updated model is evaluated; if \code{FALSE} the updated call is returned.
#' @param ... arguments to pass down.
#'
#' @importFrom stats model.matrix vcov terms predict update anova quantile weighted.mean delete.response lm lm.fit lm.wfit model.offset na.pass pchisq 
#' @importFrom lmtest coeftest waldtest waldtest.default lrtest lrtest.default
#' @importFrom car linearHypothesis
#' @import Formula sandwich

#' @rdname ivreg_Methods
#' @export
vcov.ivreg <- function(object, ...)
  object$sigma^2 * object$cov.unscaled

#' @rdname ivreg_Methods
#' @export    
bread.ivreg <- function (x, ...) 
    x$cov.unscaled * x$nobs

#' @rdname ivreg_Methods
#' @export
estfun.ivreg <- function (x, ...) 
{
    xmat <- model.matrix(x)
    if(any(alias <- is.na(coef(x)))) xmat <- xmat[, !alias, drop = FALSE]
    wts <- weights(x)
    if(is.null(wts)) wts <- 1
    res <- residuals(x)
    rval <- as.vector(res) * wts * xmat
    attr(rval, "assign") <- NULL
    attr(rval, "contrasts") <- NULL
    return(rval)
}

#' #' @rdname ivreg_Methods
#' #' @export
#' hatvalues.ivreg <- function(model, ...) {
#'   xz <- model.matrix(model, component = "projected")
#'   x  <- model.matrix(model, component = "regressors")
#'   z  <- model.matrix(model, component = "instruments")
#'   solve_qr <- function(x) chol2inv(qr.R(qr(x)))
#'   diag(x %*% solve_qr(xz) %*% t(x) %*% z %*% solve_qr(z) %*% t(z))
#' }

#' @rdname ivreg_Methods
#' @export
terms.ivreg <- function(x, component = c("regressors", "instruments", "full"), ...)
  x$terms[[match.arg(component)]]

#' @rdname ivreg_Methods
#' @export
model.matrix.ivreg <- function(object, component = c("regressors", "projected", "instruments"), ...) {
  component <- match.arg(component)
  if(!is.null(object$x)) rval <- object$x[[component]]
    else if(!is.null(object$model)) {
      X <- model.matrix(object$terms$regressors, object$model, contrasts = object$contrasts$regressors)
      Z <- if(is.null(object$terms$instruments)) NULL
        else model.matrix(object$terms$instruments, object$model, contrasts = object$contrasts$instruments)
      w <- weights(object)
      XZ <- if(is.null(Z)) {
        X
      } else if(is.null(w)) {
        lm.fit(Z, X)$fitted.values
      } else {
        lm.wfit(Z, X, w)$fitted.values
      }
      if(is.null(dim(XZ))) {
        XZ <- matrix(XZ, ncol = 1L, dimnames = list(names(XZ), colnames(X)))
	attr(XZ, "assign") <- attr(X, "assign")
      }
      rval <- switch(component,
        "regressors" = X,
	"instruments" = Z,
	"projected" = XZ)
    } else stop("not enough information in fitted model to return model.matrix")
  return(rval)
}

#' @rdname ivreg_Methods
#' @param type For \code{predict}, one of \code{"response"} (the default)  or \code{"terms"};
#' for \code{residuals}, one of \code{"response"} (the default), \code{"projected"}, \code{"regressors"},
#' \code{"working"}, \code{"deviance"}, \code{"pearson"}, or \code{"partial"}; 
#' \code{type = "working"} and \code{"response"} are equivalent, as are 
#' \code{type = "deviance"} and \code{"pearson"}.
#' 
#' @export
predict.ivreg <- function(object, newdata, type = c("response", "terms"), na.action = na.pass,  ...)
{
  type <- match.arg(type)
  if (type == "response"){
    if(missing(newdata)) fitted(object)
    else {
      mf <- model.frame(delete.response(object$terms$full), newdata,
                        na.action = na.action, xlev = object$levels)
      X <- model.matrix(delete.response(object$terms$regressors), mf,
                        contrasts = object$contrasts$regressors)
      ok <- !is.na(object$coefficients)
      drop(X[, ok, drop = FALSE] %*% object$coefficients[ok])
    } 
  } else {
    NextMethod()
  }
}

#' @rdname ivreg_Methods
#' @export
print.ivreg <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
  cat("Coefficients:\n")
  print.default(format(coef(x), digits = digits), print.gap = 2, quote = FALSE)
  cat("\n")
  invisible(x)
}

#' @rdname ivreg_Methods
#' @export
summary.ivreg <- function(object, vcov. = NULL, df = NULL, diagnostics = TRUE, ...)
{

  if (length(formula(object, component="instruments")) == 0) diagnostics <- FALSE 
      # prevent some "inherited" "lm" methods from failing
  
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

  ## covariance matrix
  if(is.null(vcov.)) 
      vc <- vcov(object)
  else {
      if(is.function(vcov.)) vc <- vcov.(object)
        else vc <- vcov.
  }
  
  ## Wald test of each coefficient
  cf <- lmtest::coeftest(object, vcov. = vc, df = df, ...)
  attr(cf, "method") <- NULL
  class(cf) <- "matrix"
  
  ## Wald test of all coefficients
  Rmat <- if(attr(object$terms$regressors, "intercept"))
    cbind(0, diag(length(na.omit(coef(object)))-1)) else diag(length(na.omit(coef(object))))
  waldtest <- car::linearHypothesis(object, Rmat, vcov. = vcov., test = ifelse(df > 0, "F", "Chisq"), singular.ok = TRUE)
  waldtest <- c(waldtest[2, "F"], waldtest[2, "Pr(>F)"], waldtest[2, "Df"], if(df > 0) waldtest[2, "Res.Df"] else NULL)
  
  ## diagnostic tests
  diag <- if(diagnostics) ivdiag(object, vcov. = vcov.) else NULL
  
  rval <- list(
    call = object$call,
    terms = object$terms,
    residuals = res,
    weights <- object$weights,
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

#' @rdname ivreg_Methods
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

#' @rdname ivreg_Methods
#' @export    
anova.ivreg <- function(object, object2, test = "F", vcov. = NULL, ...)
{
  rval <- waldtest(object, object2, test = test, vcov = vcov.)
  if(is.null(vcov.)) {
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

#' @rdname ivreg_Methods
#' @export
update.ivreg <- function (object, formula., ..., evaluate = TRUE)
{
  if(is.null(call <- getCall(object))) stop("need an object with call component")
  extras <- match.call(expand.dots = FALSE)$...
  if(!missing(formula.)) call$formula <- formula(update(Formula(formula(object)), formula.))
  if(length(extras)) {
    existing <- !is.na(match(names(extras), names(call)))
    for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
    if(any(!existing)) {
      call <- c(as.list(call), extras[!existing])
      call <- as.call(call)
    }
  }
  if(evaluate) eval(call, parent.frame())
  else call
}

ivdiag <- function(obj, vcov. = NULL) {
  ## extract data
  y <- model.response(model.frame(obj))
  x <- model.matrix(obj, component = "regressors")
  z <- model.matrix(obj, component = "instruments")
  w <- weights(obj)
  
  ## endogenous/instrument variables
  endo <- which(!(colnames(x) %in% colnames(z)))
  inst <- which(!(colnames(z) %in% colnames(x)))
  if((length(endo) <= 0L) | (length(inst) <= 0L))
    stop("no endogenous/instrument variables")

  ## return value
  rval <- matrix(NA, nrow = length(endo) + 2L, ncol = 4L)
  colnames(rval) <- c("df1", "df2", "statistic", "p-value")
  rownames(rval) <- c(if(length(endo) > 1L) paste0("Weak instruments (", colnames(x)[endo], ")") else "Weak instruments",
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

## If #Instruments = #Regressors then
##   b = (Z'X)^{-1} Z'y
## and solves the estimating equations
##   Z' (y - X beta) = 0
## For
##   cov(y) = Omega
## the following holds
##   cov(b) = (Z'X)^{-1} Z' Omega Z (X'Z)^{-1}
##   
## Generally:  
##   b = (X' P_Z X)^{-1} X' P_Z y
## with estimating equations
##   X' P_Z (y - X beta) = 0
## where P_Z is the usual projector (hat matrix wrt Z) and
##   cov(b) = (X' P_Z X)^{-1} X' P_Z Omega P_Z X (X' P_Z X)^{-1}
## Thus meat is X' P_Z Omega P_Z X and bread i (X' P_Z X)^{-1}
## 
## See
##   http://www.stata.com/support/faqs/stat/2sls.html

#' @rdname ivreg_Methods
#' @importFrom stats residuals
#' @export
residuals.ivreg <- function(object, type=c("response", "projected", "regressors", "working",
                                            "deviance", "pearson", "partial"), ...){
  type <- match.arg(type)
  w <- weights(object)
  if (is.null(w)) w <- 1
  res <- switch(type,
                working  =,
                response  = object$residuals,
                deviance =,
                pearson  = sqrt(w)*object$residuals,
                projected   = object$residuals.1,
                regressors  = object$residuals.2,
                partial  = object$residuals + predict(object, type = "terms"))
  naresid(object$na.action, res)
}

#' @rdname ivreg_Methods
#' @importFrom effects Effect
#' @param focal.predictors Focal predictors for effect plot, see \code{\link[effects]{Effect}}.
#' @export
Effect.ivreg <- function (focal.predictors, mod, ...) {
  mod$contrasts <- mod$contrasts$regressors
  NextMethod()
}

#' @rdname ivreg_Methods
#' @importFrom stats formula
#' @export
formula.ivreg <- function(x, component = c("regressors", "instruments", "complete"), ... ) {
  component <- match.arg(component)
  if (component == "complete"){
    class(x) <- "default"
    formula(x)
    
  } else {
    formula(x$terms[[component]])
  }
}

#' @rdname ivreg_Methods
#' @importFrom car Anova
#' @export
Anova.ivreg <- function(mod, test.statistic=c("F", "Chisq"), ...){
  test.statistic <- match.arg(test.statistic)
  NextMethod(test.statistic=test.statistic)
}

#' @rdname ivreg_Methods
#' @export
linearHypothesis.ivreg <- function(model, test=c("F", "Chisq"), ...){
  test <- match.arg(test)
  NextMethod(test=test)
}
