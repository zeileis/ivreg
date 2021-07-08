# Version 0.6-1

* `method` is now an explicit argument to `ivreg()` and not just passed through `...`
  to `ivreg.fit()`.
  
* More efficient computation of regression diagnostics (thanks to improvements implemented by Nikolas Kuschnig).

* Small fixes.


# Version 0.6-0

* Three-part right-hand side `formula`s are supported now to facilitate specification
  of models with many exogenous regressors. For example, if there is one exogenous
  regressor `ex` and one endogenous regressor `en` with instrument `in`, a formula
  with three parts on the right-hand side can now also be used: `y ~ ex | en | in`.
  This is equivalent to specifying: `y ~ en + ex | in + ex`.

* Robust-regression estimators are provided as an alternative to ordinary
  least squares (OLS) both in stage 1 and 2 by means of `rlm()` from package
  [MASS](https://CRAN.R-project.org/package=MASS). Specifically, in addition to
  2-stage least squares (2SLS, `method = "OLS"`, default) `ivreg()` now supports
  2-stage M-estimation (2SM, `method = "M"`) and 2-stage MM-estimation (2SMM,
  `method = "MM"`).

* Dedicated `confint()` method allowing specification of the variance-covariance
  matrix `vcov.` and degrees of freedom `df` to be used (as in the `summary()`
  method).

* Include information about which `"regressors"` are endogenous variables and
  which `"instruments"` are instruments for the endogenous variables in the
  fitted model objects from `ivreg()` and `ivreg.fit()`. Both provide elements
  `$endogenous` and `$instruments` which are named integer vectors provided
  that endogenous/instrument variables exist, and integers of length zero if
  not.
  
* Include `df.residual1` element in `ivreg` objects with the residual degrees
  of freedom from the stage-1 regression.

* Add `coef(..., component = "stage1")`, `vcov(..., component = "stage1")`, and
  `confint(..., component = "stage1")` for the estimated coefficients and
  corresponding variance-covariance matrix and confidence intervals from the
  stage-1 regression (only for the endogenous regressors). (Prompted by a request
  from Grant McDermott.)
  
* Add `residuals(..., type = "stage1")` with the residuals from the stage-1
  regression (only for the endogenous regressors).

* The `coef()`, `vcov()`, and `confint()` methods gained a `complete = TRUE` argument
  assuring that the elements pertaining to aliased coefficients are included.
  By setting `complete = FALSE` these elements are dropped.

* Include demonstration how to use `ivreg()` results in model summary tables
  and plots using the [modelsummary](https://CRAN.R-project.org/package=modelsummary)
  package.
  
* Small edits to the Diagnostics vignette.


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
