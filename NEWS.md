# Version 0.5-1

* Include information about which `"regressors"` are endogenous variables and
  which `"instruments"` are instruments for the endogenous variables in the
  fitted model objects from `ivreg()` and `ivreg.fit()`. Both provide elements
  `$endogenous` and `$instruments` which are named integer vectors provided
  that endogenous/instrument variables exist, and integers of length zero if
  not.


# Version 0.5-0

* Initial version of the `ivreg` package: An implementation of instrumental
  variables regression using two-stage least-squares (2SLS) estimation, based on
  the `ivreg()` function previously in the
  [AER](https://CRAN.R-project.org/package=AER) package. In addition to standard
  regression functionality (parameter estimation, inference, predictions, etc.)
  the package provides various regression diagnostics, including hat values,
  deletion diagnostics such as studentized residuals and Cook's distances;
  graphical diagnostics such as component-plus-residual plots and added-variable
  plots; and effect plots with partial residuals.
  
* An overview of the package, documentation, examples, and vignettes are provided
  at <https://john-d-fox.github.io/ivreg/>.
