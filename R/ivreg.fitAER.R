#' Fitting Instrumental-Variable Regressions
#' 
#' Fit instrumental-variable regression by two-stage least squares. This is
#' equivalent to direct instrumental-variables estimation when the number of
#' instruments is equal to the number of predictors.
#' 
#' \code{\link{ivreg}} is the high-level interface to the work-horse function
#' \code{ivreg.fit}. \code{ivreg.fit} is a convenience interface to \code{\link{lm.fit}} (or
#' \code{\link{lm.wfit}}) for first projecting \code{x} onto the image of
#' \code{z} and the running a regression of \code{y} onto the projected
#' \code{x}.
#' 
#' @aliases ivreg.fit
#' 
#' @param x regressor matrix.
#' @param y vector with dependent variable.
#' @param z instruments matrix.
#' @param weights an optional vector of weights to be used in the fitting
#' process.
#' @param offset an optional offset that can be used to specify an a priori
#' known component to be included during fitting.
#' @param \dots further arguments passed to \code{\link[stats:lmfit]{lm.fit}}
#' or \code{\link[stats]{lm.wfit}}, respectively.
#' @return \code{ivreg.fit} returns an unclassed list with the following
#' components: 
#' \item{coefficients}{parameter estimates, from the stage-2 regression.} 
#' \item{residuals}{vector of model residuals.} 
#' \item{residuals.1}{matrix of residuals from the stage-1 regression.}
#' \item{residuals.2}{vector of residuals from the stage-2 regression.}
#' \item{fitted.values}{vector of predicted means.}
#' \item{weights}{either the vector of weights used (if any) or \code{NULL} (if none).} 
#' \item{offset}{either the offset used (if any) or \code{NULL} (if none).} 
#' \item{estfun}{a matrix containing the empirical estimating functions.} 
#' \item{n}{number of observations.} 
#' \item{nobs}{number of observations with non-zero weights.} 
#' \item{p}{number of columns in the model matrix x of regressors.}
#' \item{q}{number of columns in the instrumental variables model matrix z}
#' \item{rank}{numeric rank of the stage-2 regression.} 
#' \item{df.residual}{residual degrees of freedom for fitted model.} 
#' \item{cov.unscaled}{unscaled covariance matrix for the coefficients.} 
#' \item{sigma}{residual standard error.}
#' \item{x}{projection of x matrix onto span of z.}
#' \item{qr}{QR decomposition for the stage-2 regression.}
#' \item{qr.1}{QR decomposition for the stage-1 regression.}
#' \item{rank.1}{numeric rank of the stage-1 regression.}
#' \item{coefficients.1}{matrix of coefficients from the stage-1 regression.}
#' @seealso \code{\link{ivreg}}, \code{\link{lm.fit}}, \code{\link{lm.wfit}}
#' @keywords regression
#' @examples
#' 
#' ## data
#' if (length(find.package("AER", quiet=TRUE)) > 0){
    #' data("CigarettesSW", package="AER")
    #' CigarettesSW$rprice <- with(CigarettesSW, price/cpi)
    #' CigarettesSW$rincome <- with(CigarettesSW, income/population/cpi)
    #' CigarettesSW$tdiff <- with(CigarettesSW, (taxs - tax)/cpi)
    #' 
    #' ## high-level interface
    #' fm <- ivreg(log(packs) ~ log(rprice) + log(rincome) | log(rincome) + tdiff + I(tax/cpi),
    #'   data = CigarettesSW, subset = year == "1995")
    #' 
    #' ## low-level interface
    #' y <- fm$y
    #' x <- model.matrix(fm, component = "regressors")
    #' z <- model.matrix(fm, component = "instruments")
    #' ivreg.fit(x, y, z)$coefficients
#' }
#' 
#' @export
ivreg.fit <- function(x, y, z, weights, offset, ...)
{
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
    # pz <- diag(NROW(x))
    # colnames(pz) <- rownames(pz) <- rownames(x)
  }
  
  ## main regression
  fit <- if(is.null(weights)) lm.fit(xz, y, offset = offset, ...)
    else lm.wfit(xz, y, weights, offset = offset, ...)
 
  ## model fit information
  ok <- which(!is.na(fit$coefficients))
  yhat <- drop(x[, ok, drop = FALSE] %*% fit$coefficients[ok])
  names(yhat) <- names(y)
  res <- y - yhat
  ucov <- chol2inv(fit$qr$qr[1:length(ok), 1:length(ok), drop = FALSE])
  colnames(ucov) <- rownames(ucov) <- names(fit$coefficients[ok])
  rss <- if(is.null(weights)) sum(res^2) else sum(weights * res^2)
  ## hat <- diag(x %*% ucov %*% t(x) %*% pz)
  ## names(hat) <- rownames(x)

  rval <- list(
    coefficients = fit$coefficients,
    residuals = res,
    residuals.1 = auxreg$residuals,
    residuals.2 = fit$residuals,
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
    sigma = sqrt(rss/fit$df.residual), ## NOTE: Stata divides by n here and uses z tests rather than t tests...
    # hatvalues = hat,
    x = xz,
    qr = fit$qr,
    qr.1 = auxreg$qr,
    rank.1 = auxreg$rank,
    coefficients.1 = coef(auxreg)
  )
  
  return(rval)
}
