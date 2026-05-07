# Summary Method for Object of Class `metadynmeta`

Summary Method for Object of Class `metadynmeta`

## Usage

``` r
# S3 method for class 'metadynmeta'
summary(
  object,
  alpha = NULL,
  ci_type = "wald",
  robust = NULL,
  nrep = 20000L,
  seed = NULL,
  ncores = NULL,
  digits = 4,
  ...
)
```

## Arguments

- object:

  an object of class `metadynmeta`.

- alpha:

  Numeric vector. Significance level \\\alpha\\. If `NULL`, the function
  will check `object` for `alpha` used in model fitting.

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

- digits:

  Integer indicating the number of decimal places to display.

- ...:

  further arguments.

## Value

Returns a matrix of estimates, standard errors, test statistics, degrees
of freedom, p-values, and confidence intervals.

## Author

Ivan Jacob Agaloos Pesigan
