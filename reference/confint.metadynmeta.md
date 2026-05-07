# Confidence Intervals for the Parameter Estimates

Confidence Intervals for the Parameter Estimates

## Usage

``` r
# S3 method for class 'metadynmeta'
confint(
  object,
  parm = NULL,
  level = 0.95,
  ci_type = "wald",
  robust = NULL,
  nrep = 20000L,
  seed = NULL,
  ncores = NULL,
  ...
)
```

## Arguments

- object:

  an object of class `metadynmeta`.

- parm:

  a specification of which parameters are to be given confidence
  intervals, either a vector of numbers or a vector of names. If
  missing, all parameters are considered.

- level:

  the confidence level required.

- ci_type:

  Character string. Valid values are `"wald"` and `"mc"`.

- robust:

  Logical. If `TRUE`, use robust (sandwich) sampling variance-covariance
  matrix. If `FALSE`, use normal theory sampling variance-covariance
  matrix. If `NULL`, the function will check `object` if robust standard
  errors are available.

- nrep:

  Positive integer. Number of replications for `ci_type = "mc"`.

- seed:

  Random seed for `ci_type = "mc"`.

- ncores:

  Positive integer. Number of cores to use for `ci_type = "mc"`.

- ...:

  further arguments.

## Value

Returns a matrix of confidence intervals.

## Author

Ivan Jacob Agaloos Pesigan
