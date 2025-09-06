# Version 0.6-6

* If a model is specified where the instruments perfectly predict all regressors,
  a warning is now issued that no endogenous variables could be detected and that
  all regressors appear to be exogenous (raised by Alex Hayes in #23).
  
* Additionally, the documentation clarifies that fitting a model with only
  a single right-hand side (i.e., treating all regressors as exogenous)
  is possible for convenience, e.g., to facilitate model comparisons. This
  now also works correctly for M and MM estimation.


# Version 0.6-5

* Better documentation for summary and inference methods for `ivreg()` objects
  with a dedicated manual page `?summary.ivreg`. Along with the `summary()`
  method this also documents the methods for `confint()`, `anova()`,
  `Anova()`, and `linearHypothesis()`. All of these take an argument `vcov.`
  so that alternative (e.g., so-called "robust") covariance matrices can be
  plugged in. The `vcov.` processing is made somewhat more convenient and
  consistent (suggested by Diogo Ferrari).


# Version 0.6-4

* Fixed bug in computing larger of stage-1 and -2 hatvalues in
  `hatvalues.ivreg()` (reported by Vasilis Syrgkanis).


# Version 0.6-3

* Enhanced `predict.ivreg()` method, which optionally provides standard errors,
  confidence intervals, and prediction intervals for predicted values.

* The `tinytable` rather than the `kableExtra` package (recently not actively maintained)
  is used now for the `modelsummary` table shown in the package vignette (contributed
  by Vincent Arel-Bundock).

* Further small improvements in the package vignettes.

* Improve non-anchored links in manual pages (prompted by CRAN).


# Version 0.6-2

* Achim Zeileis took over maintenance, both on CRAN and on GitHub. The GitHub
  source repository is now at <https://github.com/zeileis/ivreg/> with the web
  page at <https://zeileis.github.io/ivreg/>.

* Avoid partial argument matches by calling `model.matrix(..., contrasts.arg = ...)`
  rather than just `contrasts` (reported by Kevin Tappe).
  
* Make names of arguments of `influencePlot.ivreg()` and `outlierTest.ivreg()`
  consistent with the corresponding generic functions from the car package.


# Version 0.6-1

* `method` is now an explicit argument to `ivreg()` and not just passed through `...`
  to `ivreg.fit()`.
  
* More efficient computation of regression diagnostics (thanks to improvements
  implemented by Nikolas Kuschnig).

* In models without any exogenous variables (i.e., not even an exogenous `(Intercept)`)
  the `$instruments` element in the fitted model object was erroneously empty, leading
  to some incorrect subsequent computations. Also the `$endogenous` element was an
  unnamed (rather than named) vector. Both problems have been fixed now.
  (Reported by Luke Sonnet.)

* In the `summary()` method the default is now `diagnostics = NULL` (rather than
  always `TRUE`). It is now only set to `TRUE` if there are both endogenous and
  instrument variables, and `FALSE` otherwise. (Reported by Brantly Callaway.)

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
  at `https://john-d-fox.github.io/ivreg/`.
