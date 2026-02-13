# Variance-Covariance Matrix Method for an Object of Class `metadynmeta`

Variance-Covariance Matrix Method for an Object of Class `metadynmeta`

## Usage

``` r
# S3 method for class 'metadynmeta'
vcov(object, robust = NULL, ...)
```

## Arguments

- object:

  an object of class `metadynmeta`.

- robust:

  Logical. If `TRUE`, use robust (sandwich) sampling variance-covariance
  matrix. If `FALSE`, use normal theory sampling variance-covariance
  matrix. If `NULL`, the function will check `object` if robust standard
  errors are available.

- ...:

  further arguments.

## Value

Returns the sampling variance-covariance matrix of the estimated
parameters.

## Author

Ivan Jacob Agaloos Pesigan
