# Print Method for Object of Class `metadynmeta`

Print Method for Object of Class `metadynmeta`

## Usage

``` r
# S3 method for class 'metadynmeta'
print(x, alpha = NULL, robust = NULL, digits = 4, ...)
```

## Arguments

- x:

  an object of class `metadynmeta`.

- alpha:

  Numeric vector. Significance level \\\alpha\\. If `NULL`, the function
  will check `object` for `alpha` used in model fitting.

- robust:

  Logical. If `TRUE`, use robust (sandwich) sampling variance-covariance
  matrix. If `FALSE`, use normal theory sampling variance-covariance
  matrix. If `NULL`, the function will check `object` if robust standard
  errors are available.

- digits:

  Integer indicating the number of decimal places to display.

- ...:

  further arguments.

## Value

Returns a matrix of estimates, standard errors, test statistics, degrees
of freedom, p-values, and confidence intervals.

## Author

Ivan Jacob Agaloos Pesigan
