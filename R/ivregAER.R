#' Instrumental-Variable Regression
#' 
#' Fit instrumental-variable regression by two-stage least squares. This is
#' equivalent to direct instrumental-variables estimation when the number of
#' instruments is equal to the number of predictors.
#' 
#' \code{ivreg} is the high-level interface to the work-horse function
#' \code{\link{ivreg.fit}}, a set of standard methods (including \code{print},
#' \code{summary}, \code{vcov}, \code{anova}, \code{predict}, \code{residuals},
#' \code{terms}, \code{model.matrix}, \code{bread}, \code{estfun}) is available
#' and described on \code{\link{ivreg_Methods}}. For methods related to regression
#' diagnotics, see \code{\link{2SLS_Diagnostics}}.
#' 
#' Regressors and instruments for \code{ivreg} are most easily specified in a
#' formula with two parts on the right-hand side, e.g., \code{y ~ x1 + x2 | z1
#' + z2 + z3}, where \code{x1} and \code{x2} are the regressors and \code{z1},
#' \code{z2}, and \code{z3} are the instruments. Note that exogenous regressors
#' have to be included as instruments for themselves. For example, if there is
#' one exogenous regressor \code{ex} and one endogenous regressor \code{en}
#' with instrument \code{in}, the appropriate formula would be \code{y ~ ex +
#' en | ex + in}. Equivalently, this can be specified as \code{y ~ ex + en | .
#' - en + in}, i.e., by providing an update formula with a \code{.} in the
#' second part of the formula. The latter is typically more convenient, if
#' there is a large number of exogenous regressors.
#' 
#' @aliases ivreg
#' @param formula,instruments formula specification(s) of the regression
#' relationship and the instruments. Either \code{instruments} is missing and
#' \code{formula} has three parts as in \code{y ~ x1 + x2 | z1 + z2 + z3}
#' (recommended) or \code{formula} is \code{y ~ x1 + x2} and \code{instruments}
#' is a one-sided formula \code{~ z1 + z2 + z3} (only for backward
#' compatibility).
#' @param data an optional data frame containing the variables in the model.
#' By default the variables are taken from the environment of the
#' \code{formula}.
#' @param subset an optional vector specifying a subset of observations to be
#' used in fitting the model.
#' @param na.action a function that indicates what should happen when the data
#' contain \code{NA}s. The default is set by the \code{na.action} option.
#' @param weights an optional vector of weights to be used in the fitting
#' process.
#' @param offset an optional offset that can be used to specify an a priori
#' known component to be included during fitting.
#' @param contrasts an optional list. See the \code{contrasts.arg} of
#' \code{\link[stats:model.matrix]{model.matrix.default}}.
#' @param model,x,y logicals.  If \code{TRUE} the corresponding components of
#' the fit (the model frame, the model matrices , the response) are returned.
#' @param \dots further arguments passed to \code{\link{ivreg.fit}}.
#' 
#' @return \code{ivreg} returns an object of class \code{"ivreg"} that inherits from
#' class \code{"lm"}, with the following components: 
#' \item{coefficients}{parameter estimates.}
#' \item{residuals}{a vector of residuals.} 
#' \item{fitted.values}{a vector of predicted means.} 
#' \item{weights}{either the vector of weights used (if any) or \code{NULL} (if none).} 
#' \item{offset}{either the offset used (if any) or \code{NULL} (if none).} 
#' \item{n}{number of observations.} 
#' \item{nobs}{number of observations with non-zero weights.} 
#' \item{rank}{the numeric rank of the fitted linear model.} 
#' \item{df.residual}{residual degrees of freedom for fitted model.} 
#' \item{cov.unscaled}{unscaled covariance matrix for the coefficients.} 
#' \item{sigma}{residual standard error.} 
#' \item{call}{the original function call.} 
#' \item{formula}{the model formula.}
#' \item{na.action}{function applied to missing values in the model fit.}
#' \item{terms}{a list with elements \code{"regressors"} and \code{"instruments"} 
#' containing the terms objects for the respective components.} 
#' \item{levels}{levels of the categorical regressors.} 
#' \item{contrasts}{the contrasts used for categorical regressors.} 
#' \item{model}{the full model frame (if \code{model = TRUE}).} 
#' \item{y}{the response vector (if \code{y = TRUE}).} 
#' \item{x}{a list with elements \code{"regressors"}, \code{"instruments"}, \code{"projected"},
#' containing the model matrices from the respective components (if \code{x = TRUE}). 
#' \code{"projected"} is the matrix of regressors projected on the image of the instruments.}
#' @seealso \code{\link{ivreg.fit}}, \code{\link{2SLS_Diagnostics}}, \code{\link{ivreg_Methods}},
#' \code{\link[stats]{lm}}, \code{\link[stats:lmfit]{lm.fit}}
#' @references Greene, W. H. (1993) \emph{Econometric Analysis}, 2nd ed.,
#' Macmillan.
#' @keywords regression
#' @examples
#' 
#' ## data
#' if (length(find.package("AER", quiet=TRUE)) > 0){
    #' data("CigarettesSW", package = "AER")
    #' CigarettesSW$rprice <- with(CigarettesSW, price/cpi)
    #' CigarettesSW$rincome <- with(CigarettesSW, income/population/cpi)
    #' CigarettesSW$tdiff <- with(CigarettesSW, (taxs - tax)/cpi)
    #' 
    #' ## model 
    #' fm <- ivreg(log(packs) ~ log(rprice) + log(rincome) | log(rincome) + tdiff + I(tax/cpi),
    #'   data = CigarettesSW, subset = year == "1995")
    #' summary(fm)
    #' summary(fm, vcov = sandwich::sandwich, df = Inf, tests = TRUE)
    #' 
    #' ## ANOVA
    #' fm2 <- ivreg(log(packs) ~ log(rprice) | tdiff, data = CigarettesSW, subset = year == "1995")
    #' anova(fm, fm2)
#' }
#' 
#' @export
ivreg <- function(formula, instruments, data, subset, na.action, weights, offset,
  contrasts = NULL, model = TRUE, y = TRUE, x = FALSE, ...)
{
  ## set up model.frame() call  
  cl <- match.call()
  if(missing(data)) data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action", "weights", "offset"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  
  ## handle instruments for backward compatibility
  if(!missing(instruments)) {
    formula <- Formula::as.Formula(formula, instruments)
    cl$instruments <- NULL
    cl$formula <- formula(formula)
  } else {
    formula <- Formula::as.Formula(formula)
  }
  stopifnot(length(formula)[1] == 1L, length(formula)[2] %in% 1:2)
  
  ## try to handle dots in formula
  has_dot <- function(formula) inherits(try(terms(formula), silent = TRUE), "try-error")
  if(has_dot(formula)) {
    f1 <- formula(formula, rhs = 1)
    f2 <- formula(formula, lhs = 0, rhs = 2)
    if(!has_dot(f1) & has_dot(f2)) formula <- Formula::as.Formula(f1,
      update(formula(formula, lhs = 0, rhs = 1), f2))
  }
  
  ## call model.frame()
  mf$formula <- formula
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  
  ## extract response, terms, model matrices
  Y <- model.response(mf, "numeric")
  mt <- terms(formula, data = data)
  mtX <- terms(formula, data = data, rhs = 1)
  X <- model.matrix(mtX, mf, contrasts)
  if(length(formula)[2] < 2L) {
    mtZ <- NULL
    Z <- NULL
  } else {
    mtZ <- delete.response(terms(formula, data = data, rhs = 2))
    Z <- model.matrix(mtZ, mf, contrasts)
  }

  ## weights and offset
  weights <- model.weights(mf)
  offset <- model.offset(mf)
  if(is.null(offset)) offset <- 0
  if(length(offset) == 1) offset <- rep(offset, NROW(Y))
  offset <- as.vector(offset)

  ## call default interface
  rval <- ivreg.fit(X, Y, Z, weights, offset, ...)

  ## enhance information stored in fitted model object
  rval$call <- cl
  rval$formula <- formula(formula)
  rval$terms <- list(regressors = mtX, instruments = mtZ, full = mt)
  rval$na.action <- attr(mf, "na.action")
  rval$levels <- .getXlevels(mt, mf)
  rval$contrasts <- list(regressors = attr(X, "contrasts"), instruments = attr(Z, "contrasts"))
  if(model) rval$model <- mf
  if(y) rval$y <- Y
  if(x) rval$x <- list(regressors = X, instruments = Z, projected = rval$x)
    else rval$x <- NULL
      
  class(rval) <- c("ivreg", "lm")
  return(rval)
}
