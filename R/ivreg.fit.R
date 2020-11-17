#' Fitting Instrumental-Variable Regressions by 2SLS
#' 
#' Fit instrumental-variable regression by two-stage least squares (2SLS). This is
#' equivalent to direct instrumental-variables estimation when the number of
#' instruments is equal to the number of predictors.
#' 
#' \code{\link{ivreg}} is the high-level interface to the work-horse function
#' \code{ivreg.fit}. \code{ivreg.fit} is essentially a convenience interface to
#' \code{\link[stats:lmfit]{lm.fit}} (or \code{\link[stats:lmfit]{lm.wfit}})
#' for first projecting \code{x} onto the image of
#' \code{z}, then running a regression of \code{y} on the projected
#' \code{x}, and computing the residual standard deviation.
#' 
#' @aliases ivreg.fit
#' 
#' @importFrom MASS rlm
#' @importFrom stats mad
#' 
#' @param x regressor matrix.
#' @param y vector for the response variable.
#' @param z instruments matrix.
#' @param weights an optional vector of weights to be used in the fitting
#' process.
#' @param offset an optional offset that can be used to specify an a priori
#' known component to be included during fitting.
#' @param method the method used to fit the stage 1 and 2 regression: 
#' \code{"OLS"} for traditional 2SLS regression (the default), 
#' \code{"M"} for M-estimation, or \code{"MM"} for MM-estimation, with the
#' latter two robust-regression methods implemented via the \code{\link[MASS]{rlm}} 
#' function in the \pkg{MASS} package.
#' @param rlm.args a list of optional arguments to be passed to the \code{\link[MASS]{rlm}} 
#' function in the \pkg{MASS} package if robust regression is used for the stage 1 and 2 regressions.
#' @param \dots further arguments passed to \code{\link[stats:lmfit]{lm.fit}}
#' or \code{\link[stats:lmfit]{lm.wfit}}, respectively.
#' @return \code{ivreg.fit} returns an unclassed list with the following
#' components: 
#' \item{coefficients}{parameter estimates, from the stage-2 regression.} 
#' \item{residuals}{vector of model residuals.} 
#' \item{residuals1}{matrix of residuals from the stage-1 regression.}
#' \item{residuals2}{vector of residuals from the stage-2 regression.}
#' \item{fitted.values}{vector of predicted means for the response.}
#' \item{weights}{either the vector of weights used (if any) or \code{NULL} (if none).} 
#' \item{offset}{either the offset used (if any) or \code{NULL} (if none).} 
#' \item{estfun}{a matrix containing the empirical estimating functions.} 
#' \item{n}{number of observations.} 
#' \item{nobs}{number of observations with non-zero weights.} 
#' \item{p}{number of columns in the model matrix x of regressors.}
#' \item{q}{number of columns in the instrumental variables model matrix z}
#' \item{rank}{numeric rank of the model matrix for the stage-2 regression.} 
#' \item{df.residual}{residual degrees of freedom for fitted model.} 
#' \item{cov.unscaled}{unscaled covariance matrix for the coefficients.} 
#' \item{sigma}{residual standard error; when method is \code{"M"} or \code{"MM"}, this
#' is based on the MAD of the residuals (around 0) --- see \code{\link[stats]{mad}}.}
#' \item{x}{projection of x matrix onto span of z.}
#' \item{qr}{QR decomposition for the stage-2 regression.}
#' \item{qr1}{QR decomposition for the stage-1 regression.}
#' \item{rank1}{numeric rank of the model matrix for the stage-1 regression.}
#' \item{coefficients1}{matrix of coefficients from the stage-1 regression.}
#' \item{df.residual1}{residual degrees of freedom for the stage-1 regression.} 
#' \item{exogenous}{columns of the \code{"regressors"} matrix that are exogenous.}
#' \item{endogenous}{columns of the \code{"regressors"} matrix that are endogenous.}
#' \item{instruments}{columns of the \code{"instruments"} matrix that are
#' instruments for the endogenous variables.}
#' \item{method}{the method used for the stage 1 and 2 regressions, one of \code{"OLS"},
#' \code{"M"}, or \code{"MM"}.}
#' \item{rweights}{a matrix of robustness weights with columns for each of the stage-1
#' regressions and for the stage-2 regression (in the last column) if the fitting method is 
#' \code{"M"} or \code{"MM"}, \code{NULL} if the fitting method is \code{"OLS"}.}
#' \item{hatvalues}{a matrix of hatvalues. For \code{method = "OLS"}, the matrix consists of two
#' columns, for each of the stage-1 and stage-2 regression; for \code{method = "M"} or \code{"MM"},
#' there is one column for \emph{each} stage=1 regression and for the stage-2 regression. }
#' @seealso \code{\link{ivreg}}, \code{\link[stats:lmfit]{lm.fit}}, 
#' \code{\link[stats:lmfit]{lm.wfit}}, \code{\link[MASS]{rlm}}, \code{\link[stats]{mad}}
#' @keywords regression
#' @examples
#' ## data
#' data("CigaretteDemand", package = "ivreg")
#' 
#' ## high-level interface
#' m <- ivreg(log(packs) ~ log(rprice) + log(rincome) | salestax + log(rincome),
#'   data = CigaretteDemand)
#' 
#' ## low-level interface
#' y <- m$y
#' x <- model.matrix(m, component = "regressors")
#' z <- model.matrix(m, component = "instruments")
#' ivreg.fit(x, y, z)$coefficients
#' 
#' @export
ivreg.fit <- function(x, y, z, weights, offset, method = c("OLS", "M", "MM"),
                      rlm.args=list(), ...)
{
  
  method <- match.arg(method, c("OLS", "M", "MM"))
  
  ## model dimensions
  n <- NROW(y)
  p <- ncol(x)
  
  ## defaults
  if(missing(z)) z <- NULL
  if(missing(weights)) weights <- NULL
  if(missing(offset)) offset <- rep(0, n)
  
  ## sanity checks
  stopifnot(n == nrow(x))
  if(!is.null(z)) stopifnot(n == nrow(z))
  if(!is.null(weights)) stopifnot(n == NROW(weights))
  stopifnot(n == NROW(offset))
  
  ## project regressors x on image of instruments z
  if(!is.null(z)) {
    if(ncol(z) < ncol(x)) warning("more regressors than instruments")
    auxreg <- if(is.null(weights)) lm.fit(z, x, ...) else lm.wfit(z, x, weights, ...)
    xz <- as.matrix(auxreg$fitted.values)
    # pz <- z %*% chol2inv(auxreg$qr$qr) %*% t(z)
    colnames(xz) <- colnames(x)
  } else {
    auxreg <- NULL
    xz <- x
  }
  
  if (method == "OLS"){
    hats <- if (!is.null(auxreg)) cbind(lm.influence(auxreg, do.coef = FALSE)$hat, 0)
    else cbind(NA, matrix(0, nrow=n, ncol=1))
  }
  
  ## infer endogenous variables in x and instruments in z
  
  exog <- structure(seq_along(colnames(x)), .Names = colnames(x))
  if(!is.null(auxreg)) {
    endo <- which(colMeans(as.matrix(auxreg$residuals^2)) > sqrt(.Machine$double.eps))
    inst <- which(rowMeans(as.matrix(coef(auxreg)^2)[, -endo, drop = FALSE]) < sqrt(.Machine$double.eps))
    exog <- exog[-endo]
  } else {
    endo <- inst <- integer()
  }
  
  # robust regression for stage 1
  if (method != "OLS" && length(endo) > 0){
    rlm.args$x <- z
    rlm.args$method <- method
    if (!is.null(weights)) rlm.args$weights <- weights
    residuals <- matrix(0, n, p)
    hats <- rwts <- matrix(0, n, length(endo) + 1)
    coef <- coef(auxreg)
    j <- 0
    for (en in endo){
      j <- j + 1
      rlm.args$y <- x[, en]
      xz[, en] <- fitted(st1 <- do.call(rlm, rlm.args))
      residuals[ , en] <- residuals(st1)
      coef[, en] <- coef(st1)
      rwts[, j] <- st1$w
      hats[, j] <- hatvalues(st1)
    } 
    auxreg$residuals <- residuals
    auxreg$coefficients <- coef
  }
  
  ## main regression
  fit <- if (method == "OLS") {
    if(is.null(weights)) lm.fit(xz, y, offset = offset, ...)
    else lm.wfit(xz, y, weights, offset = offset, ...)
  } else {
    # robust regression for stage 2
    rlm.args$y <- y - offset
    rlm.args$x <- xz
    do.call(rlm, rlm.args)
  }
  if (method != "OLS") {
    rwts[, ncol(rwts)] <- fit$w
    hats[, ncol(hats)] <- hatvalues(fit)
    rownames(hats) <- rownames(rwts) <- names(y)
    colnames(hats) <- colnames(rwts) <- c(names(endo), "stage_2")
    fit$df.residual <- n - length(na.omit(coef(fit)))
  } else {
    hats[, 2] <- lm.influence(fit, do.coef = FALSE)$hat
    colnames(hats) <- c("stage_1", "stage_2")
  }
  
  ## model fit information
  ok <- which(!is.na(fit$coefficients))
  yhat <- drop(x[, ok, drop = FALSE] %*% fit$coefficients[ok]) + offset
  names(yhat) <- names(y)
  res <- y - yhat
  ucov <- chol2inv(fit$qr$qr[1:length(ok), 1:length(ok), drop = FALSE])
  colnames(ucov) <- rownames(ucov) <- names(fit$coefficients[ok])
  sigma <- if (method == "OLS") {
    rss <- if(is.null(weights)) sum(res^2) else sum(weights * res^2)
    sqrt(rss/fit$df.residual) ## NOTE: Stata divides by n here and uses z tests rather than t tests...
  } else {
    if (is.null(weights)) mad(res, 0)  else mad(sqrt(weights) * res, 0)
  }

  rval <- list(
    coefficients = fit$coefficients,
    residuals = res,
    residuals1 = auxreg$residuals,
    residuals2 = fit$residuals,
    fitted.values = yhat,
    weights = weights,
    offset = if(identical(offset, rep(0, n))) NULL else offset,
    n = n,
    nobs = if(is.null(weights)) n else sum(weights > 0),
    p = p,
    q = ncol(z),
    rank = fit$rank,
    df.residual = fit$df.residual,
    cov.unscaled = ucov,
    sigma = sigma,
    # hatvalues = hat,
    x = xz,
    qr = fit$qr,
    qr1 = auxreg$qr,
    rank1 = auxreg$rank,
    coefficients1 = coef(auxreg),
    df.residual1 = auxreg$df.residual,
    exogenous = exog,
    endogenous = endo,
    instruments = inst,
    method = method,
    rweights = if (method == "OLS") NULL else rwts,
    hatvalues = hats
  )
  
  return(rval)
}
