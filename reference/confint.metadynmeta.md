# Confidence Intervals for the Parameter Estimates

Confidence Intervals for the Parameter Estimates

## Usage

``` r
# S3 method for class 'metadynmeta'
confint(object, parm = NULL, level = 0.95, lb = TRUE, robust = FALSE, ...)
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

- lb:

  Logical. If `TRUE`, returns profile likelihood-based confidence
  intervals. If `FALSE`, returns Wald confidence intervals.

- robust:

  Logical. If `TRUE`, use robust (sandwich) sampling variance-covariance
  matrix. If `FALSE`, use normal theory sampling variance-covariance
  matrix.

- ...:

  further arguments.

## Value

Returns a matrix of confidence intervals.

## Author

Ivan Jacob Agaloos Pesigan
