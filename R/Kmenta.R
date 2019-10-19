#' Partly Artifical Data on the U. S. Economy
#'
#' @docType data
#'
#' @usage Kmenta
#'
#' @description These are partly contrived data from Kmenta (1986), constructed
#' to illustrate estimation of a simultaneous-equation econometric model. The data
#' are an annual time-series for the U.S. economy from 1922 to 1941. The values of the
#' exogenous variables D, and F, and A are real, while those of the endogenous
#' variables Q and P are simulated according to the linear simultaneous equation model 
#' fit in the examples.
#'
#' @format A data frame with 20 rows and 5 columns.
#' \describe{
#'   \item{Q}{food consumption per capita.}
#'   \item{P}{ratio of food prices to general consumer prices.}
#'   \item{D}{disposible income in constant dollars.}
#'   \item{F}{ratio of preceding year's prices received by farmers to general consumer prices.}
#'   \item{A}{time in years.}
#'   }
#'
#' @source Kmenta, J. (1986) \emph{Elements of Econometrics, Second Edition}, Macmillan.
#' @seealso \code{\link{ivreg}}.
#'
#' @examples 
#' deq <- ivreg(Q ~ P + D | D + F + A, data=Kmenta) # demand equation
#' summary(deq, diagnostics=TRUE)
#' seq <- ivreg(Q ~ P + F + A | D + F + A, data=Kmenta) # supply equation
#' summary(seq, diagnostics=TRUE)
"Kmenta"
