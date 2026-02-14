# Meta-Regression with Distal Outcome

## Dynamics Description

The *Stable Persistence with Baseline-Moderated Autoregression* process
represents a bivariate dynamic system in which two latent psychological
constructs (e.g., positive and negative affect) exhibit within-construct
persistence over time through autoregressive dynamics. The transition
matrix is diagonal, implying that each construct evolves according to
its own self-regulatory process and there are no cross-lagged
(reciprocal) influences between the constructs in the systematic
dynamics.

Between-person heterogeneity in persistence is captured via a
time-invariant baseline covariate $`x_i \in \left\{ 0, 1 \right\}`$ with
one value per individual. The individual-specific transition matrix is
modeled as
``` math
\begin{equation}
  \mathrm{vec} \left( \boldsymbol{\beta}_i \right) = \mathrm{vec} \left( \boldsymbol{\beta} \left( x_i \right) \right) = \mathrm{vec} \left( \boldsymbol{\beta}_0 \right) + \mathrm{vec} \left( \boldsymbol{\beta}_1 \right) x_i,
\end{equation}
```
so that individuals with $`x_i = 0`$ follows $`\boldsymbol{\beta}_0`$,
whereas individuals with $`x_i = 1`$ follow
$`\boldsymbol{\beta}_0 + \boldsymbol{\beta}_1`$. Under the current
specification, $`\boldsymbol{\beta}_1`$ increases the diagonal
autoregressive parameters, implying that the $`x_i = 1`$ group exhibits
greater persistence—deviations from an individual’s equilibrium decay
more slowly—while remaining dynamically stable.

Between-person heterogeneity in baseline levels is also captured via the
same time-invariant covariate $`x_i \in \left\{ 0, 1 \right\}`$ through
the person-specific set-point vector
``` math
\begin{equation}
  \boldsymbol{\mu}_i = \boldsymbol{\mu} \left( x_i \right) = \boldsymbol{\mu}_{0} + \boldsymbol{\mu}_{1} x_i, 
\end{equation}
```
so that individuals with $`x_i = 0`$ follow $`\boldsymbol{\mu}_{0}`$,
whereas individuals with $`x_i = 1`$ follow
$`\boldsymbol{\mu}_{0} + \boldsymbol{\mu}_{1}`$.

In addition to explaining between-person differences in the DT-VAR
parameters via the baseline covariate $`x_i`$, we also include a distal
outcome $`z_i`$ measured once per individual. Conceptually, $`z_i`$ can
be interpreted as a later outcome (e.g., symptom severity, functioning,
or substance use) whose between-person differences are partly explained
by (a) each person’s estimated set-point and dynamic parameters and (b)
$`x_i`$.

The process noise covariance allows for small disturbances that may be
correlated across constructs, permitting coordinated innovations even
though the lagged dynamics are decoupled. Measurement errors are assumed
to be minimal and symmetric across indicators.

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
estimates the set-point vector $`\boldsymbol{\mu}_{i}`$ is given by
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

### Distal outcome model

Let $`z_i`$ denote a distal outcome observed once per individual $`i`$.
We model $`z_i`$ as a linear function of (a) the person-specific DT-VAR
parameters and (b) the baseline covariate $`x_i`$:
``` math
\begin{equation}
  z_i = \kappa + \boldsymbol{\phi}^{\prime} \mathbf{w}_i + \omega x_i + \xi_i, \quad \mathrm{with} \quad \xi_i \sim \mathcal{N}(0, \psi),
\end{equation}
```
where $`\kappa`$ is an intercept, $`\boldsymbol{\phi}`$ is a vector of
regression coefficients, $`\omega`$ is the direct effect of $`x_i`$ on
$`z_i`$, and $`\psi`$ is the residual variance of the distal outcome.
The predictor vector $`\mathbf{w}_i`$ stacks the person-specific DT-VAR
coefficients and intercepts used for prediction. In this vignette we use
``` math
\begin{equation}
  \mathbf{w}_i = \left[ \mathrm{vec}(\boldsymbol{\beta}_i)^{\prime}, \; \boldsymbol{\nu}_i^{\prime} \right]^{\prime},
\end{equation}
```
so $`\mathbf{w}_i`$ has length $`p^2 + p`$ (here $`4 + 2 = 6`$).

## Data Generation

### Notation

Let $`t = 100`$ be the number of time points and $`n = 1000`$ be the
number of individuals. We simulate a total of time \$= 1.01 &times;
10\<sup\>4\</sup\>\$ points per individual, discarding the first
$`10<sup>4</sup>`$ as burn-in. The analysis uses the final $`100`$
measurement occasions.

Let the initial condition $`\boldsymbol{\eta}_{0}`$ be given by
``` math
\begin{equation}
  \boldsymbol{\eta}_{0} \sim \mathcal{N} \left( \boldsymbol{\mu}_{\boldsymbol{\eta} \mid 0}, \boldsymbol{\Sigma}_{\boldsymbol{\eta} \mid 0} \right) .
\end{equation}
```
$`\boldsymbol{\mu}_{\boldsymbol{\eta} \mid 0}`$ and
$`\boldsymbol{\Sigma}_{\boldsymbol{\eta} \mid 0}`$ are functions of
$`\boldsymbol{\alpha}`$ and $`\boldsymbol{\beta}`$.

Let the intercept vector $`\boldsymbol{\alpha}`$ when $`X = 0`$ be
``` math
\begin{equation}
  \left(
    \begin{array}{c}
      0.5 \\
      -0.5 \\
    \end{array}
  \right) .
\end{equation}
```
Let the intercept vector $`\boldsymbol{\alpha}`$ when $`X = 1`$ be
``` math
\begin{equation}
\left(
\begin{array}{c}
  0.75 \\
  -0.75 \\
\end{array}
\right) .
\end{equation}
```
Let the covariance matrix be
``` math
\begin{equation}
  \left(
    \begin{array}{cc}
      0.1 & -0.05 \\
      -0.05 & 0.1 \\
    \end{array}
  \right) .
\end{equation}
```

Let the transition matrix $`\boldsymbol{\beta}`$ when $`X = 0`$ be
``` math
\begin{equation}
  \left(
    \begin{array}{cc}
      0.5 & 0 \\
      0 & 0.5 \\
    \end{array}
  \right) .
\end{equation}
```
Let the transition matrix $`\boldsymbol{\beta}`$ when $`X = 1`$ be
``` math
\begin{equation}
\left(
\begin{array}{cc}
  0.75 & 0 \\
  0 & 0.75 \\
\end{array}
\right) .
\end{equation}
```
Let the covariance matrix be and covariance matrix
``` math
\begin{equation}
  \left(
    \begin{array}{cccc}
      0.02 & 0.01 & 0 & 0 \\
      0.01 & 0.015 & 0 & 0 \\
      0 & 0 & 0.01 & 0.005 \\
      0 & 0 & 0.005 & 0.015 \\
    \end{array}
  \right) .
\end{equation}
```

The `SimAlphaN` and `SimBetaNCovariate` functions from the
`simStateSpace` package generate random intercept vectors and transition
matrices from the multivariate normal distribution. Note that the
`SimBetaNCovariate` function generates transition matrices that are
weakly stationary with an option to set lower and upper bounds. The
person-specific set-point vector $`\boldsymbol{\mu}_{i}`$ was derived
from the generated $`\boldsymbol{\alpha}_{i}`$ and
$`\boldsymbol{\beta}_{i}`$.

Let the dynamic process noise $`\boldsymbol{\Psi}`$ be given by
``` math
\begin{equation}
  \boldsymbol{\Psi}
  =
  \left(
    \begin{array}{cc}
      0.2 & -0.05 \\
      -0.05 & 0.18 \\
    \end{array}
  \right) .
\end{equation}
```

### R Function Arguments

``` r

n
#> [1] 1000
time
#> [1] 10100
burnin
#> [1] 10000
# first mu0 in the list of length n
mu0[[1]]
#> [1]  1.4001203 -0.4648662
# first sigma0 in the list of length n
sigma0[[1]]
#>            [,1]       [,2]
#> [1,]  0.2518659 -0.0466736
#> [2,] -0.0466736  0.2516021
# first sigma0_l in the list of length n
sigma0_l[[1]] # sigma0_l <- t(chol(sigma0))
#>             [,1]      [,2]
#> [1,]  0.50186241 0.0000000
#> [2,] -0.09300079 0.4929026
# first alpha in the list of length n
alpha[[1]]
#> [1]  0.7602558 -0.4388544
# first beta in the list of length n
beta[[1]]
#>           [,1]        [,2]
#> [1,] 0.4479308 -0.02733575
#> [2,] 0.1605012  0.53936547
# first psi in the list of length n
psi[[1]]
#>       [,1]  [,2]
#> [1,]  0.20 -0.05
#> [2,] -0.05  0.18
psi_l[[1]] # psi_l <- t(chol(psi))
#>            [,1]      [,2]
#> [1,]  0.4472136 0.0000000
#> [2,] -0.1118034 0.4092676
# first mu (set-point) in the list of length n
mu[[1]]
#> [1]  1.4001203 -0.4648662
# distal outcome parameters
kappa_z
#> [1] 10
phi_z
#>      [,1] [,2] [,3] [,4] [,5] [,6]
#> [1,]    5    5    5    5    5    5
omega_z
#> [1] 0.15
psi_z
#> [1] 0.45
```

### Visualizing the Dynamics Without Process Noise and Measurement Error (n = 5 with Different Initial Condition)

#### $`X = 0`$

![](fig-vignettes-distal-no-error-0-1.png)![](fig-vignettes-distal-no-error-0-2.png)

#### $`X = 1`$

![](fig-vignettes-distal-no-error-1-1.png)![](fig-vignettes-distal-no-error-1-2.png)

### Using the `SimSSMIVary` Function from the `simStateSpace` Package to Simulate Data

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
data <- as.data.frame(sim, burnin = burnin)
head(data)
#>   id time       y1          y2
#> 1  1    0 2.348797 -0.36133130
#> 2  1    1 2.203791 -0.06306478
#> 3  1    2 1.763803 -0.11966576
#> 4  1    3 1.608552 -0.88481869
#> 5  1    4 1.116157 -0.36118972
#> 6  1    5 2.254971 -0.37188289
plot(sim, burnin = burnin)
```

![](fig-vignettes-distal-error-1.png)![](fig-vignettes-distal-error-2.png)

## Model Fitting

``` r

library(OpenMx)
library(fitVARMxID)
```

The `FitVARMxID` function fits a VAR model on each individual $`i`$.

### LDL’-parameterized covariance matrices

Covariances such as `psi` and `theta` are estimated using the LDL’
decomposition of a positive definite covariance matrix. The
decomposition expresses a covariance matrix $`\Sigma`$ as\
``` math
\begin{equation}
  \boldsymbol{\Sigma} = \left( \mathbf{L} + \mathbf{I} \right) \mathrm{diag} \left( \mathrm{Softplus} \left( \mathbf{d}_{uc} \right) \right) \left( \mathbf{L} + \mathbf{I} \right)^{\prime},
\end{equation}
```
where: - $`\mathbf{L}`$ is a strictly lower-triangular matrix of free
parameters (`l_mat_strict`),\
- $`\mathbf{I}`$ is the identity matrix,\
- $`\mathbf{d}_{uc}`$ is an unconstrained vector,\
-
$`\mathrm{Softplus} \left(\mathbf{d}_{uc} \right) = \log \left(1 + \exp \left( \mathbf{d}_{uc} \right) \right)`$
ensures strictly positive diagonal entries.

The
[`LDL()`](https://github.com/jeksterslab/fitVARMxID/reference/LDL.html)
function extracts this decomposition from a positive definite covariance
matrix. It returns:\
- `d_uc`: unconstrained diagonal parameters, equal to
`InvSoftplus(d_vec)`,\
- `d_vec`: diagonal entries, equal to `Softplus(d_uc)`,\
- `l_mat_strict`: the strictly lower-triangular factor.

``` r

sigma <- matrix(
  data = c(1.0, 0.5, 0.5, 1.0),
  nrow = 2,
  ncol = 2
)
ldl_sigma <- LDL(sigma)
d_uc <- ldl_sigma$d_uc
l_mat_strict <- ldl_sigma$l_mat_strict
I <- diag(2)
sigma_reconstructed <- (l_mat_strict + I) %*% diag(log1p(exp(d_uc)), 2) %*% t(l_mat_strict + I)
sigma_reconstructed
#>      [,1] [,2]
#> [1,]  1.0  0.5
#> [2,]  0.5  1.0
```

### `FitVARMxID`

``` r

fit <- FitVARMxID(
  data = data,
  observed = c("y1", "y2"),
  id = "id",
  center = TRUE,
  tries_explore = 1000,
  tries_local = 100,
  max_attempts = 100,
  ncores = parallel::detectCores()
)
```

#### Parameter estimates

``` r

head(summary(fit))
#>                               mu[1,1]    mu[2,1] beta[1,1]   beta[2,1]
#> FitVARMxID_DTVAR_ID3.Rds  1.153348406 -0.5873101 0.6935418  0.17246835
#> FitVARMxID_DTVAR_ID6.Rds  0.001610978 -0.4100247 0.6153076 -0.11361758
#> FitVARMxID_DTVAR_ID7.Rds  0.874685699 -1.3716365 0.4763569  0.12629573
#> FitVARMxID_DTVAR_ID9.Rds  0.120502387 -1.0458112 0.2191061 -0.14986574
#> FitVARMxID_DTVAR_ID10.Rds 1.150417875 -0.4283811 0.6651179  0.08401120
#> FitVARMxID_DTVAR_ID11.Rds 3.065915464 -1.5910143 0.8072554 -0.04304906
#>                             beta[1,2] beta[2,2]  psi[1,1]    psi[2,1]  psi[2,2]
#> FitVARMxID_DTVAR_ID3.Rds  -0.18577947 0.4557184 0.1649165 -0.02709214 0.1662770
#> FitVARMxID_DTVAR_ID6.Rds   0.13996673 0.4095755 0.2172841 -0.07375769 0.1743473
#> FitVARMxID_DTVAR_ID7.Rds  -0.11021823 0.5943148 0.1707274 -0.03526842 0.1709325
#> FitVARMxID_DTVAR_ID9.Rds   0.01075629 0.5836428 0.1388499 -0.04100114 0.1745615
#> FitVARMxID_DTVAR_ID10.Rds  0.09156199 0.3895123 0.2505575 -0.07955761 0.2248892
#> FitVARMxID_DTVAR_ID11.Rds -0.12771113 0.6032254 0.2243289 -0.04499742 0.2176939
```

#### Proportion of converged cases

``` r

converged(
  fit,
  prop = TRUE
)
#> [1] 0.468
```

## Mixed-Effects Meta-Analysis with Distal Outcome of Person-Specific Dynamics and Means

    #> Warning in mapply(FUN = function(x, y) {: longer argument not a multiple of
    #> length of shorter

We synthesize the person-specific estimates to recover population-level
effects, between-person variability, and systematic covariate-related
differences. We fit a mixed-effects meta-analytic model in which each
individual’s estimate is weighted by its within-person sampling
uncertainty, random effects capture residual heterogeneity across
individuals, and fixed effects quantify the association between
covariate $`X`$ and the person-specific dynamic parameters as well as
effects on a distal outcome.

``` r

library(metaDyn)
random <- MetaVARMx(
  fit,
  x = x,
  z = z,
  effects = TRUE,
  set_point = TRUE,
  robust_v = FALSE,
  robust = TRUE,
  lb = TRUE,
  ncores = parallel::detectCores()
)
#> Error in `[.data.frame`:
#> ! undefined columns selected
```

``` r

summary(random)
#> Error:
#> ! object 'random' not found
```

### Normal Theory Confidence Intervals

``` r

confint(random, level = 0.95, lb = FALSE)
#> Error:
#> ! object 'random' not found
```

``` r

confint(random, level = 0.99, lb = FALSE)
#> Error:
#> ! object 'random' not found
```

### Robust Confidence Intervals

``` r

confint(random, level = 0.95, lb = FALSE, robust = TRUE)
#> Error:
#> ! object 'random' not found
```

``` r

confint(random, level = 0.99, lb = FALSE, robust = TRUE)
#> Error:
#> ! object 'random' not found
```

### Profile-Likelihood Confidence Intervals

``` r

confint(random, level = 0.95, lb = TRUE)
#> Error:
#> ! object 'random' not found
```

``` r

confint(random, level = 0.99, lb = TRUE)
#> Error:
#> ! object 'random' not found
```

- The fixed part of the random-effects model gives pooled means for the
  person-specific parameters (e.g.,
  $`\mathbb{E}[\mathrm{vec}(\boldsymbol{\beta}_i)]`$ and
  $`\mathbb{E}[\boldsymbol{\nu}_i]`$) and fixed effects for how $`X`$
  shifts these parameters.
- The random part yields between-person covariances
  ($`\boldsymbol{\tau}^2`$) quantifying residual heterogeneity in
  dynamics and intercepts beyond what is explained by $`X`$.
- The distal-outcome sub-model yields regression effects linking the the
  estimated person-specific parameters (via $`\boldsymbol{\phi}`$) to
  distal outcome $`z_i`$, a direct effect of $`X`$ on $`z_i`$ (via
  $`\omega`$), and the residual variance $`\psi`$.

``` r

means <- extract(random, what = "alpha")
#> Error:
#> ! object 'random' not found
means
#> Error:
#> ! object 'means' not found
covariances <- extract(random, what = "tau_sqr")
#> Error:
#> ! object 'random' not found
covariances
#> Error:
#> ! object 'covariances' not found
```

Finally, we compare the meta-analytic population estimates to the known
generating values.

``` r

pop_mean
#> [1]  2.225402237 -2.150926350  0.613120730 -0.008650001 -0.004652319
#> [6]  0.617315267
pop_cov
#>             [,1]        [,2]         [,3]          [,4]          [,5]
#> [1,]  8.14262645  0.86433780  0.267330532  0.0622197590 -0.0946237191
#> [2,]  0.86433780  9.23200499 -0.042748641  0.1354952004 -0.0377962010
#> [3,]  0.26733053 -0.04274864  0.030048825  0.0078692525 -0.0010470459
#> [4,]  0.06221976  0.13549520  0.007869252  0.0139002330 -0.0006657224
#> [5,] -0.09462372 -0.03779620 -0.001047046 -0.0006657224  0.0096798022
#> [6,]  0.08444775 -0.24227808  0.012450141 -0.0011007600  0.0037814441
#>              [,6]
#> [1,]  0.084447745
#> [2,] -0.242278079
#> [3,]  0.012450141
#> [4,] -0.001100760
#> [5,]  0.003781444
#> [6,]  0.026568250
```

## Summary

This vignette demonstrates a two-stage hierarchical estimation approach
for dynamic systems: 1. Individual-level DT-VAR estimation. 2.
Population-level meta-analysis of person-specific dynamics and means. 3.
Estimation and interpretation of covariate effects, where $`X`$ predicts
systematic between-person differences in dynamics and baseline levels.
4. A distal outcome model linking $`z_i`$ to person-specific
dynamic/mean features and to $`X`$.

## References

Hunter, M. D. (2017). State space modeling in an open source, modular,
structural equation modeling environment. *Structural Equation Modeling:
A Multidisciplinary Journal*, *25*(2), 307–324.
<https://doi.org/10.1080/10705511.2017.1369354>

Neale, M. C., Hunter, M. D., Pritikin, J. N., Zahery, M., Brick, T. R.,
Kirkpatrick, R. M., Estabrook, R., Bates, T. C., Maes, H. H., & Boker,
S. M. (2015). OpenMx 2.0: Extended structural equation and statistical
modeling. *Psychometrika*, *81*(2), 535–549.
<https://doi.org/10.1007/s11336-014-9435-8>

R Core Team. (2024). *R: A language and environment for statistical
computing*. R Foundation for Statistical Computing.
<https://www.R-project.org/>
