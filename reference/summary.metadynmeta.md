# Summary Method for Object of Class `metadynmeta`

Summary Method for Object of Class `metadynmeta`

## Usage

``` r
# S3 method for class 'metadynmeta'
summary(object, alpha = 0.05, lb = FALSE, robust = FALSE, digits = 4, ...)
```

## Arguments

- object:

  an object of class `metadynmeta`.

- alpha:

  Numeric vector. Significance level \\\alpha\\.

- lb:

  Logical. If `TRUE`, returns profile likelihood-based confidence
  intervals. If `FALSE`, returns Wald confidence intervals.

- robust:

  Logical. If `TRUE`, use robust (sandwich) sampling variance-covariance
  matrix. If `FALSE`, use normal theory sampling variance-covariance
  matrix.

- digits:

  Integer indicating the number of decimal places to display.

- ...:

  further arguments.

## Value

Returns a matrix of estimates, standard errors, test statistics, degrees
of freedom, p-values, and confidence intervals.

## Author

Ivan Jacob Agaloos Pesigan
