# Meta-Analysis of Discrete-Time Vector Autoregressive Model Estimates

## Dynamics Description

The *Stable Reciprocal Regulation* process represents a bivariate
dynamic system in which two latent psychological constructs—such as
negative and positive affect—mutually influence each other over time.
Each construct shows moderate self-regulation (autoregressive effects)
and mild opposing cross-effects, reflecting an equilibrium-seeking
mechanism characteristic of emotional balance.

Individuals vary in their self-regulatory tendencies and in the strength
of these antagonistic couplings. At the population level, the transition
matrix indicates that increases in one construct are followed by slight
decreases in the other, producing a stable, damped oscillatory pattern
around individual equilibrium points. The process noise covariance
allows for small correlated disturbances, while measurement errors are
assumed to be minimal and symmetric across indicators.

This dynamic pattern captures a psychologically plausible process of
*reciprocal inhibition*—where short-term fluctuations in one system
component (e.g., negative affect) are naturally counteracted by
adjustments in its counterpart (e.g., positive affect), leading to
emotional homeostasis over time.

## Model

The measurement model is given by
``` math
\begin{equation}
  \mathbf{y}_{i, t} = \boldsymbol{\eta}_{i, t}
\end{equation}
```
where $`\mathbf{y}_{i, t}`$ and $`\boldsymbol{\eta}_{i, t}`$ are random
variables.

The dynamic structure is given by
``` math
\begin{equation}
  \boldsymbol{\eta}_{i, t} = \boldsymbol{\alpha}_{i} + \boldsymbol{\beta}_{i} \boldsymbol{\eta}_{i, t - 1} + \boldsymbol{\zeta}_{i, t}, \quad \mathrm{with} \quad \boldsymbol{\zeta}_{i, t} \sim \mathcal{N}   \left( \mathbf{0}, \boldsymbol{\Psi} \right)
\end{equation}
```
where $`\boldsymbol{\eta}_{i, t}`$, $`\boldsymbol{\eta}_{i, t - 1}`$,
and $`\boldsymbol{\zeta}_{i, t}`$ are random variables, and
$`\boldsymbol{\beta}_{i}`$, and $`\boldsymbol{\Psi}`$ are model
parameters. Here, $`\boldsymbol{\eta}_{i, t}`$ is a vector of latent
variables at time $`t`$ and individual $`i`$,
$`\boldsymbol{\eta}_{i, t - 1}`$ represents a vector of latent variables
at time $`t - 1`$ and individual $`i`$, and
$`\boldsymbol{\zeta}_{i, t}`$ represents a vector of dynamic noise at
time $`t`$ and individual $`i`$. $`\boldsymbol{\beta}_{i}`$ is a matrix
of autoregression and cross regression coefficients for individual
$`i`$, and $`\boldsymbol{\Psi}`$ the covariance matrix of
$`\boldsymbol{\zeta}_{i, t}`$ that is invariant across all individuals.
In this model, $`\boldsymbol{\Psi}`$ is a symmetric matrix.

### Alternative Parameterization

An alternative parameterization of the dynamic structure that directly
estimates the set point vector $`\boldsymbol{\mu}_{i}`$ is given by
``` math
\begin{equation}
  \boldsymbol{\eta}_{i, t} = \boldsymbol{\mu}_{i} + \boldsymbol{\beta}_{i} \left( \boldsymbol{\eta}_{i, t - 1} - \boldsymbol{\mu}_{i} \right) + \boldsymbol{\zeta}_{i, t} .
\end{equation}
```

Algebraic manipulation of the equation results in the following
``` math
\begin{equation}
  \boldsymbol{\eta}_{i, t} = \boldsymbol{\mu}_{i} - \boldsymbol{\beta}_{i} \boldsymbol{\mu}_{i} + \boldsymbol{\beta}_{i} \boldsymbol{\eta}_{i, t - 1} + \boldsymbol{\zeta}_{i, t} ,
\end{equation}
```
where we can see that the intercept vector $`\boldsymbol{\alpha}_{i}`$
is implied by
$`\boldsymbol{\mu}_{i} - \boldsymbol{\beta}_{i} \boldsymbol{\mu}_{i}`$.

## Data Generation

### Notation

Let $`t = 1000`$ be the number of time points and $`n = 100`$ be the
number of individuals.

Let the initial condition $`\boldsymbol{\eta}_{0}`$ be given by
``` math
\begin{equation}
  \boldsymbol{\eta}_{0} \sim \mathcal{N} \left( \boldsymbol{\mu}_{\boldsymbol{\eta} \mid 0}, \boldsymbol{\Sigma}_{\boldsymbol{\eta} \mid 0} \right) .
\end{equation}
```
$`\boldsymbol{\mu}_{\boldsymbol{\eta} \mid 0}`$ and
$`\boldsymbol{\Sigma}_{\boldsymbol{\eta} \mid 0}`$ are functions of
$`\boldsymbol{\alpha}`$ and $`\boldsymbol{\beta}`$.

Let the set point vector $`\boldsymbol{\mu}`$ be normally distributed
with the following means
``` math
\begin{equation}
  \left(
    \begin{array}{c}
      2.87 \\
      2.04 \\
    \end{array}
  \right)
\end{equation}
```
and covariance matrix
``` math
\begin{equation}
  \left(
    \begin{array}{cc}
      1.2 & 0.46 \\
      0.46 & 1.1 \\
    \end{array}
  \right) .
\end{equation}
```

Let the transition matrix $`\boldsymbol{\beta}`$ be normally distributed
with the following means
``` math
\begin{equation}
  \left(
    \begin{array}{cc}
      0.28 & -0.048 \\
      -0.035 & 0.26 \\
    \end{array}
  \right)
\end{equation}
```
and covariance matrix
``` math
\begin{equation}
  \left(
    \begin{array}{cccc}
      0.017 & 0.005 & 0.001 & 0.008 \\
      0.005 & 0.008 & 0.001 & 0.006 \\
      0.001 & 0.001 & 0.001 & 0.002 \\
      0.008 & 0.006 & 0.002 & 0.026 \\
    \end{array}
  \right) .
\end{equation}
```

The `SimMuN` and `SimBetaN` functions from the `simStateSpace` package
generate random set point vectors and transition matrices from the
multivariate normal distribution. Note that the `SimBetaN` function
generates transition matrices that are weakly stationary with an option
to set lower and upper bounds.

Let the dynamic process noise $`\boldsymbol{\Psi}`$ be given by
``` math
\begin{equation}
  \boldsymbol{\Psi}
  =
  \left(
    \begin{array}{cc}
      1.3 & 0.57 \\
      0.57 & 1.56 \\
    \end{array}
  \right) .
\end{equation}
```

### R Function Arguments

``` r

n
#> [1] 100
time
#> [1] 1000
# first mu0 in the list of length n
mu0[[1]]
#> [1] 2.287769 2.606098
# first sigma0 in the list of length n
sigma0[[1]]
#>           [,1]      [,2]
#> [1,] 1.4403549 0.5646226
#> [2,] 0.5646226 1.6885649
# first sigma0_l in the list of length n
sigma0_l[[1]] # sigma0_l <- t(chol(sigma0))
#>           [,1]     [,2]
#> [1,] 1.2001479 0.000000
#> [2,] 0.4704608 1.211293
# first alpha in the list of length n
alpha[[1]]
#> [1] 1.678968 2.010503
# first beta in the list of length n
beta[[1]]
#>             [,1]        [,2]
#> [1,]  0.32874206 -0.05498069
#> [2,] -0.07362662  0.29317209
# first psi in the list of length n
psi[[1]]
#>      [,1] [,2]
#> [1,] 1.30 0.57
#> [2,] 0.57 1.56
psi_l[[1]] # psi_l <- t(chol(psi))
#>           [,1]     [,2]
#> [1,] 1.1401754 0.000000
#> [2,] 0.4999231 1.144586
# first mu (set point) in the list of length n
mu[[1]]
#> [1] 2.287769 2.606098
```

### Visualizing the Dynamics Without Process Noise and Measurement Error (n = 5 with Different Initial Condition)

![](fig-vignettes-meta-dt-var-no-error-1.png)![](fig-vignettes-meta-dt-var-no-error-2.png)

### Using the `SimSSMVARIVary` Function from the `simStateSpace` Package to Simulate Data

``` r

library(simStateSpace)
sim <- SimSSMVARIVary(
  n = n,
  time = time,
  mu0 = mu0,
  sigma0_l = sigma0_l,
  alpha = alpha,
  beta = beta,
  psi_l = psi_l
)
data <- as.data.frame(sim)
head(data)
#>   id time        y1       y2
#> 1  1    0 2.2494783 1.541912
#> 2  1    1 2.8646254 2.346091
#> 3  1    2 0.8218736 3.348071
#> 4  1    3 2.4900816 2.282659
#> 5  1    4 2.7074414 3.944741
#> 6  1    5 2.9755157 2.667683
plot(sim)
```

![](fig-vignettes-meta-dt-var-error-1.png)![](fig-vignettes-meta-dt-var-error-2.png)

## Stage 1: Person-Specific VAR Model

``` r

library(OpenMx)
library(fitVARMxID)
```

The `FitVARMxID` function fits a VAR model on each individual $`i`$.

> **Note:** Consider using the argument `ncores` to use multiple cores
> for parallel processing.

``` r

fit <- FitVARMxID(
  data = data,
  observed = c("y1", "y2"),
  id = "id",
  center = TRUE,
  ncores = parallel::detectCores()
)
```

``` r

summary(
  fit,
  means = TRUE
)
#> Call:
#> FitVARMxID(data = data, observed = c("y1", "y2"), id = "id", 
#>     center = TRUE, ncores = parallel::detectCores())
#> 
#> Convergence:
#> 100.0%
#> 
#> Means of the estimated paramaters per individual.
#>   mu[1,1]   mu[2,1] beta[1,1] beta[2,1] beta[1,2] beta[2,2]  psi[1,1]  psi[2,1] 
#>    2.9733    2.1842    0.2711   -0.0659   -0.0543    0.2394    1.2932    0.5612 
#>  psi[2,2] 
#>    1.5390
```

## Stage 2: Random-Effects Meta-Analysis of Person-Specific Set Points and Dynamics

We synthesize the person-specific estimates to recover population-level
effects and their between-person variability. We use a random-effects
model so the pooled mean reflects both within-person estimation
uncertainty and between-person heterogeneity.

All available parameters are meta-analyzed by default. Setting
`cov_dyn = FALSE`, meta-analyzes only the set points and transition
matrix. Setting `tau_sqr_l_free`, such that covariances between `mu` and
`beta` are constained to zero, simplifies the random effects.

``` r

tau_sqr_l_free <- matrix(
  data = c(
    FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
    TRUE, FALSE, FALSE, FALSE, FALSE, FALSE,
    FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
    FALSE, FALSE, TRUE, FALSE, FALSE, FALSE,
    FALSE, FALSE, TRUE, TRUE, FALSE, FALSE,
    FALSE, FALSE, TRUE, TRUE, TRUE, FALSE
  ),
  byrow = TRUE,
  nrow = 6,
  ncol = 6
)
# FALSE values in tau_sqr_l_free will be treated as fixed values.
# TRUE values in tau_sqr_l_free will be treated as starting values.
tau_sqr_l_values <- matrix(
  data = 0,
  nrow = 6,
  ncol = 6
)
```

> **Note:** Consider using the argument `ncores` to use multiple cores
> for parallel processing.

``` r

library(metaDyn)
metavar <- MetaVARMx(
  object = fit,
  tau_sqr_l_free = tau_sqr_l_free,
  tau_sqr_l_values = tau_sqr_l_values,
  cov_dyn = FALSE,
  ncores = parallel::detectCores()
)
```

``` r

summary(metavar)
#> Call:
#> MetaVARMx(object = fit, tau_sqr_l_free = tau_sqr_l_free, tau_sqr_l_values = tau_sqr_l_values, 
#>     cov_dyn = FALSE, ncores = parallel::detectCores())
#> 
#> Status code:
#> 0
#> 
#> CI type:
#> "normal"
#> 
#>                  est     se         z      p    2.5%   97.5%
#> alpha[1,1]    2.9733 0.1115   26.6781 0.0000  2.7549  3.1918
#> alpha[2,1]    2.1844 0.1052   20.7594 0.0000  1.9781  2.3906
#> alpha[3,1]    0.2720 0.0122   22.3157 0.0000  0.2481  0.2959
#> alpha[4,1]   -0.0659 0.0106   -6.2428 0.0000 -0.0866 -0.0452
#> alpha[5,1]   -0.0545 0.0047  -11.6738 0.0000 -0.0636 -0.0453
#> alpha[6,1]    0.2401 0.0150   16.0294 0.0000  0.2107  0.2694
#> tau_sqr[1,1]  1.2396 0.1757    7.0559 0.0000  0.8952  1.5839
#> tau_sqr[2,1]  0.5389 0.1291    4.1756 0.0000  0.2859  0.7918
#> tau_sqr[2,2]  1.1040 0.1566    7.0514 0.0000  0.7971  1.4108
#> tau_sqr[3,3]  0.0138 0.0021    6.6038 0.0000  0.0097  0.0179
#> tau_sqr[4,3]  0.0021 0.0013    1.5812 0.1138 -0.0005  0.0046
#> tau_sqr[5,3]  0.0013 0.0006    2.3450 0.0190  0.0002  0.0025
#> tau_sqr[6,3]  0.0037 0.0019    2.0081 0.0446  0.0001  0.0074
#> tau_sqr[4,4]  0.0099 0.0016    6.2656 0.0000  0.0068  0.0130
#> tau_sqr[5,4]  0.0007 0.0005    1.4928 0.1355 -0.0002  0.0017
#> tau_sqr[6,4]  0.0053 0.0017    3.2207 0.0013  0.0021  0.0086
#> tau_sqr[5,5]  0.0013 0.0003    4.2040 0.0000  0.0007  0.0019
#> tau_sqr[6,5]  0.0023 0.0007    3.1304 0.0017  0.0008  0.0037
#> tau_sqr[6,6]  0.0214 0.0032    6.7778 0.0000  0.0152  0.0276
#> i_sqr[1,1]    0.9982 0.0003 3860.4969 0.0000  0.9977  0.9987
#> i_sqr[2,1]    0.9978 0.0003 3177.0151 0.0000  0.9972  0.9984
#> i_sqr[3,1]    0.9293 0.0100   93.3768 0.0000  0.9098  0.9488
#> i_sqr[4,1]    0.8879 0.0159   55.9008 0.0000  0.8568  0.9190
#> i_sqr[5,1]    0.5975 0.0572   10.4451 0.0000  0.4854  0.7096
#> i_sqr[6,1]    0.9537 0.0065  146.4238 0.0000  0.9409  0.9665
```

### Normal Theory Confidence Intervals

``` r

confint(metavar, level = 0.95)
#>                      2.5 %       97.5 %
#> alpha[1,1]    2.7548840130  3.191767598
#> alpha[2,1]    1.9781275976  2.390593226
#> alpha[3,1]    0.2481422065  0.295927269
#> alpha[4,1]   -0.0866287758 -0.045230607
#> alpha[5,1]   -0.0636118759 -0.045322487
#> alpha[6,1]    0.2107036915  0.269408544
#> tau_sqr[1,1]  0.8952340791  1.583878123
#> tau_sqr[2,1]  0.2859461519  0.791842021
#> tau_sqr[2,2]  0.7971048330  1.410797663
#> tau_sqr[3,3]  0.0097047686  0.017896770
#> tau_sqr[4,3] -0.0004944567  0.004622886
#> tau_sqr[5,3]  0.0002204809  0.002465308
#> tau_sqr[6,3]  0.0000894128  0.007368261
#> tau_sqr[4,4]  0.0067964394  0.012984098
#> tau_sqr[5,4] -0.0002307961  0.001705876
#> tau_sqr[6,4]  0.0020859982  0.008571934
#> tau_sqr[5,5]  0.0006881582  0.001890227
#> tau_sqr[6,5]  0.0008491006  0.003692931
#> tau_sqr[6,6]  0.0151913841  0.027551416
#> i_sqr[1,1]    0.9976654900  0.998679029
#> i_sqr[2,1]    0.9971648879  0.998395989
#> i_sqr[3,1]    0.9097638232  0.948774245
#> i_sqr[4,1]    0.8567713523  0.919033674
#> i_sqr[5,1]    0.4853556950  0.709578194
#> i_sqr[6,1]    0.9409422848  0.966474114
```

``` r

confint(metavar, level = 0.99)
#>                      0.5 %       99.5 %
#> alpha[1,1]    2.6862446272  3.260406983
#> alpha[2,1]    1.9133245508  2.455396273
#> alpha[3,1]    0.2406346291  0.303434846
#> alpha[4,1]   -0.0931328995 -0.038726484
#> alpha[5,1]   -0.0664853471 -0.042449016
#> alpha[6,1]    0.2014804908  0.278631744
#> tau_sqr[1,1]  0.7870402584  1.692071944
#> tau_sqr[2,1]  0.2064641504  0.871324022
#> tau_sqr[2,2]  0.7006867015  1.507215794
#> tau_sqr[3,3]  0.0084177119  0.019183827
#> tau_sqr[4,3] -0.0012984495  0.005426879
#> tau_sqr[5,3] -0.0001322070  0.002817996
#> tau_sqr[6,3] -0.0010541771  0.008511851
#> tau_sqr[4,4]  0.0058242878  0.013956250
#> tau_sqr[5,4] -0.0005350693  0.002010149
#> tau_sqr[6,4]  0.0010669838  0.009590948
#> tau_sqr[5,5]  0.0004992996  0.002079085
#> tau_sqr[6,5]  0.0004023024  0.004139730
#> tau_sqr[6,6]  0.0132494824  0.029493317
#> i_sqr[1,1]    0.9975062516  0.998838267
#> i_sqr[2,1]    0.9969714679  0.998589409
#> i_sqr[3,1]    0.9036348417  0.954903227
#> i_sqr[4,1]    0.8469892326  0.928815794
#> i_sqr[5,1]    0.4501277873  0.744806102
#> i_sqr[6,1]    0.9369309437  0.970485455
```

### Robust Confidence Intervals

``` r

confint(metavar, level = 0.95, robust = TRUE)
#>                      2.5 %       97.5 %
#> alpha[1,1]    2.7537803903  3.192871220
#> alpha[2,1]    1.9770923273  2.391628496
#> alpha[3,1]    0.2480754094  0.295994066
#> alpha[4,1]   -0.0867440323 -0.045115351
#> alpha[5,1]   -0.0636813829 -0.045252980
#> alpha[6,1]    0.2105511338  0.269561101
#> tau_sqr[1,1]  0.9013971675  1.577715034
#> tau_sqr[2,1]  0.2794426273  0.798345545
#> tau_sqr[2,2]  0.8278232911  1.380079204
#> tau_sqr[3,3]  0.0091150653  0.018486474
#> tau_sqr[4,3] -0.0002357525  0.004364182
#> tau_sqr[5,3]  0.0003160392  0.002369750
#> tau_sqr[6,3] -0.0005935359  0.008051209
#> tau_sqr[4,4]  0.0072591176  0.012521420
#> tau_sqr[5,4] -0.0001960170  0.001671097
#> tau_sqr[6,4]  0.0023413351  0.008316597
#> tau_sqr[5,5]  0.0007204290  0.001857956
#> tau_sqr[6,5]  0.0011175412  0.003424491
#> tau_sqr[6,6]  0.0134984518  0.029244348
#> i_sqr[1,1]    0.9976745608  0.998669958
#> i_sqr[2,1]    0.9972265109  0.998334366
#> i_sqr[3,1]    0.9069556478  0.951582421
#> i_sqr[4,1]    0.8614269820  0.914378044
#> i_sqr[5,1]    0.4913753405  0.703558549
#> i_sqr[6,1]    0.9374453527  0.969971046
```

``` r

confint(metavar, level = 0.99, robust = TRUE)
#>                      0.5 %       99.5 %
#> alpha[1,1]    2.684794e+00  3.261857389
#> alpha[2,1]    1.911964e+00  2.456756848
#> alpha[3,1]    2.405468e-01  0.303522632
#> alpha[4,1]   -9.328437e-02 -0.038575011
#> alpha[5,1]   -6.657669e-02 -0.042357668
#> alpha[6,1]    2.012800e-01  0.278832239
#> tau_sqr[1,1]  7.951399e-01  1.683972272
#> tau_sqr[2,1]  1.979171e-01  0.879871102
#> tau_sqr[2,2]  7.410576e-01  1.466844897
#> tau_sqr[3,3]  7.642710e-03  0.019958829
#> tau_sqr[4,3] -9.584545e-04  0.005086884
#> tau_sqr[5,3] -6.622135e-06  0.002692411
#> tau_sqr[6,3] -1.951724e-03  0.009409397
#> tau_sqr[4,4]  6.432350e-03  0.013348187
#> tau_sqr[5,4] -4.893618e-04  0.001964442
#> tau_sqr[6,4]  1.402553e-03  0.009255379
#> tau_sqr[5,5]  5.417105e-04  0.002036674
#> tau_sqr[6,5]  7.550931e-04  0.003786939
#> tau_sqr[6,6]  1.102459e-02  0.031718208
#> i_sqr[1,1]    9.975182e-01  0.998826346
#> i_sqr[2,1]    9.970525e-01  0.998508423
#> i_sqr[3,1]    8.999443e-01  0.958593795
#> i_sqr[4,1]    8.531078e-01  0.922697259
#> i_sqr[5,1]    4.580389e-01  0.736894947
#> i_sqr[6,1]    9.323352e-01  0.975081203
```

- The fixed part of the random-effects model gives pooled means
  $`\boldsymbol{\alpha} = \mathbb{E} \left[ \mathrm{Vec} \left( \boldsymbol{\mu}, \boldsymbol{\beta} \right)  \right]`$.
- The random part yields between-person covariances
  ($`\boldsymbol{\tau}^{2}`$) quantifying heterogeneity in set point
  ($`\boldsymbol{\mu}`$) and dynamics ($`\boldsymbol{\beta}`$) across
  individuals.

``` r

means <- extract(metavar, what = "alpha")
means
#>          alpha
#> y1  2.97332581
#> y2  2.18436041
#> y3  0.27203474
#> y4 -0.06592969
#> y5 -0.05446718
#> y6  0.24005612
covariances <- extract(metavar, what = "tau_sqr")
covariances
#>           y1        y2          y3          y4          y5          y6
#> y1 1.2395561 0.5388941 0.000000000 0.000000000 0.000000000 0.000000000
#> y2 0.5388941 1.1039512 0.000000000 0.000000000 0.000000000 0.000000000
#> y3 0.0000000 0.0000000 0.013800769 0.002064215 0.001342895 0.003728837
#> y4 0.0000000 0.0000000 0.002064215 0.009890269 0.000737540 0.005328966
#> y5 0.0000000 0.0000000 0.001342895 0.000737540 0.001289192 0.002271016
#> y6 0.0000000 0.0000000 0.003728837 0.005328966 0.002271016 0.021371400
```

## Summary

This vignette demonstrates a two-stage hierarchical estimation approach
for dynamic systems: 1. individual-level DT-VAR estimation, and\
2. population-level meta-analysis of person-specific set points and
dynamics.

## References

Cheung, M. W.-L. (2015). *Meta-analysis: A structural equation modeling
approach*. Wiley. <https://doi.org/10.1002/9781118957813>

Hunter, M. D. (2017). State space modeling in an open source, modular,
structural equation modeling environment. *Structural Equation Modeling:
A Multidisciplinary Journal*, *25*(2), 307–324.
<https://doi.org/10.1080/10705511.2017.1369354>

Neale, M. C., Hunter, M. D., Pritikin, J. N., Zahery, M., Brick, T. R.,
Kirkpatrick, R. M., Estabrook, R., Bates, T. C., Maes, H. H., & Boker,
S. M. (2015). OpenMx 2.0: Extended structural equation and statistical
modeling. *Psychometrika*, *81*(2), 535–549.
<https://doi.org/10.1007/s11336-014-9435-8>

R Core Team. (2026). *R: A language and environment for statistical
computing*. R Foundation for Statistical Computing.
<https://www.R-project.org/>
