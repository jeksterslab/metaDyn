# Meta-Regression

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

The process noise covariance allows for small disturbances that may be
correlated across constructs, permitting coordinated innovations even
though the lagged dynamics are decoupled. Measurement errors are assumed
to be minimal and symmetric across indicators.

## Model

The measurement model is given by
``` math
\begin{equation}
  \mathbf{y}_{i, t} = \boldsymbol{\Lambda} \boldsymbol{\eta}_{i, t} + \boldsymbol{\varepsilon}_{i, t}, \quad \mathrm{with} \quad \boldsymbol{\varepsilon}_{i, t} \sim \mathcal{N} \left( \mathbf{0}, \boldsymbol{\Theta} \right)
\end{equation}
```
where $`\mathbf{y}_{i, t}`$, $`\boldsymbol{\eta}_{i, t}`$, and
$`\boldsymbol{\varepsilon}_{i, t}`$ are random variables and
$`\boldsymbol{\Lambda}`$, and $`\boldsymbol{\Theta}`$ are model
parameters. $`\mathbf{y}_{i, t}`$ represents a vector of observed random
variables, $`\boldsymbol{\eta}_{i, t}`$ a vector of latent random
variables, and $`\boldsymbol{\varepsilon}_{i, t}`$ a vector of random
measurement errors, at time $`t`$ and individual $`i`$.
$`\boldsymbol{\Lambda}`$ denotes a matrix of factor loadings, and
$`\boldsymbol{\Theta}`$ the covariance matrix of
$`\boldsymbol{\varepsilon}`$ that is invariant across individuals. In
this model, $`\boldsymbol{\Lambda}`$ is an identity matrix and
$`\boldsymbol{\Theta}`$ is a symmetric matrix.

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

## Data Generation

### Notation

Let $`t = 1000`$ be the number of time points and $`n = 1000`$ be the
number of individuals. We simulate a total of time $`= 11000`$ points
per individual, discarding the first $`10000`$ as burn-in. The analysis
uses the final $`1000`$ measurement occasions.

Let the factor loadings matrix $`\boldsymbol{\Lambda}`$ be given by
``` math
\begin{equation}
  \boldsymbol{\Lambda}
  =
  \left(
    \begin{array}{cc}
      1 & 0 \\
      0 & 1 \\
    \end{array}
  \right) .
\end{equation}
```

Let the measurement error covariance matrix $`\boldsymbol{\Theta}`$ be
given by
``` math
\begin{equation}
  \boldsymbol{\Theta}
  =
  \left(
  \begin{array}{cc}
    0.5 & 0 \\
    0 & 0.5 \\
  \end{array}
  \right) .
\end{equation}
```

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
#> [1] 11000
burnin
#> [1] 10000
# first mu0 in the list of length n
mu0[[1]]
#> [1]  0.7520376 -0.5127468
# first sigma0 in the list of length n
sigma0[[1]]
#>            [,1]       [,2]
#> [1,]  0.2787159 -0.1023574
#> [2,] -0.1023574  0.2538853
# first sigma0_l in the list of length n
sigma0_l[[1]] # sigma0_l <- t(chol(sigma0))
#>            [,1]      [,2]
#> [1,]  0.5279355 0.0000000
#> [2,] -0.1938823 0.4650752
# first alpha in the list of length n
alpha[[1]]
#> [1]  0.3319244 -0.1910529
# first beta in the list of length n
beta[[1]]
#>             [,1]       [,2]
#> [1,]  0.47646139 -0.1205201
#> [2,] -0.08920883  0.4965522
# first psi in the list of length n
psi[[1]]
#>       [,1]  [,2]
#> [1,]  0.20 -0.05
#> [2,] -0.05  0.18
psi_l[[1]] # psi_l <- t(chol(psi))
#>            [,1]      [,2]
#> [1,]  0.4472136 0.0000000
#> [2,] -0.1118034 0.4092676
nu
#> [[1]]
#> [1] 0 0
lambda
#> [[1]]
#>      [,1] [,2]
#> [1,]    1    0
#> [2,]    0    1
theta
#> [[1]]
#>      [,1] [,2]
#> [1,]  0.5  0.0
#> [2,]  0.0  0.5
theta_l # theta_l <- t(chol(theta))
#> [[1]]
#>           [,1]      [,2]
#> [1,] 0.7071068 0.0000000
#> [2,] 0.0000000 0.7071068
# first mu_eta (set-point) in the list of length n
mu_eta[[1]]
#> [1]  0.7520376 -0.5127468
```

### Visualizing the Dynamics Without Process Noise and Measurement Error (n = 5 with Different Initial Condition)

#### $`X = 0`$

![](fig-vignettes-covariate-no-error-0-1.png)![](fig-vignettes-covariate-no-error-0-2.png)

#### $`X = 1`$

![](fig-vignettes-covariate-no-error-1-1.png)![](fig-vignettes-covariate-no-error-1-2.png)

### Using the `SimSSMIVary` Function from the `simStateSpace` Package to Simulate Data

``` r

library(simStateSpace)
sim <- SimSSMIVary(
  n = n,
  time = time,
  mu0 = mu0,
  sigma0_l = sigma0_l,
  alpha = alpha,
  beta = beta,
  psi_l = psi_l,
  nu = nu,
  lambda = lambda,
  theta_l = theta_l
)
data <- as.data.frame(sim, burnin = burnin)
head(data)
#>   id time         y1         y2
#> 1  1    0 -0.2039857 -0.5271213
#> 2  1    1  2.3635709 -0.8225477
#> 3  1    2  1.0336287 -1.3055546
#> 4  1    3  0.3748494 -0.1758500
#> 5  1    4 -0.6424827 -1.8117723
#> 6  1    5  1.0244405 -0.1976583
plot(sim, burnin = burnin)
```

![](fig-vignettes-covariate-error-1.png)![](fig-vignettes-covariate-error-2.png)

## Model Fitting

``` r

library(OpenMx)
library(fitDTVARMxID)
```

The `FitDTVARMxID` function fits a DT-VAR model on each individual
$`i`$. To set up the estimation, we first provide **starting values**
for each parameter matrix.

### Set-Point (`mu_eta`)

The set-point vector $`\boldsymbol{\mu}`$ is initialized with starting
values.

``` r

mu_eta_values <- mu_eta
```

### Autoregressive Parameters (`beta`)

We initialize the autoregressive coefficient matrix
$`\boldsymbol{\beta}`$ with the true values used in simulation.

``` r

beta_values <- beta
```

### LDL’-parameterized covariance matrices

Covariances such as `psi` and `theta` are estimated using the LDL’
decomposition of a positive definite covariance matrix. The
decomposition expresses a covariance matrix $`\Sigma`$ as  
``` math
\begin{equation}
  \boldsymbol{\Sigma} = \left( \mathbf{L} + \mathbf{I} \right) \mathrm{diag} \left( \mathrm{Softplus} \left( \mathbf{d}_{uc} \right) \right) \left( \mathbf{L} + \mathbf{I} \right)^{\prime},
\end{equation}
```
where: - $`\mathbf{L}`$ is a strictly lower-triangular matrix of free
parameters (`l_mat_strict`),  
- $`\mathbf{I}`$ is the identity matrix,  
- $`\mathbf{d}_{uc}`$ is an unconstrained vector,  
-
$`\mathrm{Softplus} \left(\mathbf{d}_{uc} \right) = \log \left(1 + \exp \left( \mathbf{d}_{uc} \right) \right)`$
ensures strictly positive diagonal entries.

The
[`LDL()`](https://github.com/jeksterslab/fitDTVARMxID/reference/LDL.html)
function extracts this decomposition from a positive definite covariance
matrix. It returns:  
- `d_uc`: unconstrained diagonal parameters, equal to
`InvSoftplus(d_vec)`,  
- `d_vec`: diagonal entries, equal to `Softplus(d_uc)`,  
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

#### Process Noise Covariance Matrix (`psi`)

Starting values for the process noise covariance matrix
$`\boldsymbol{\Psi}`$ are given below, with corresponding LDL’
parameters.

``` r

psi_values <- psi[[1]]
ldl_psi_values <- LDL(psi_values)
psi_d_values <- ldl_psi_values$d_uc
psi_l_values <- ldl_psi_values$l_mat_strict
```

``` r

psi_d_values
#> [1] -1.507772 -1.701853
```

``` r

psi_l_values
#>       [,1] [,2]
#> [1,]  0.00    0
#> [2,] -0.25    0
```

#### Measurement Error Covariance Matrix (`theta`)

Starting values for the measurement error covariance matrix
$`\boldsymbol{\Theta}`$ are given below, with corresponding LDL’
parameters.

``` r

theta_values <- theta[[1]]
ldl_theta_values <- LDL(theta_values)
theta_d_values <- ldl_theta_values$d_uc
theta_l_values <- ldl_theta_values$l_mat_strict
```

``` r

theta_d_values
#> [1] -0.4327521 -0.4327521
```

``` r

theta_l_values
#>      [,1] [,2]
#> [1,]    0    0
#> [2,]    0    0
```

### Initial mean vector (`mu_0`) and covariance matrix (`sigma_0`)

The initial mean vector $`\boldsymbol{\mu_0}`$ and covariance matrix
$`\boldsymbol{\Sigma_0}`$ are fixed using `mu0` and `sigma0`.

``` r

mu0_values <- mu0
```

``` r

sigma0_values <- lapply(
  X = sigma0,
  FUN = LDL
)
sigma0_d_values <- lapply(
  X = sigma0_values,
  FUN = function(i) {
    i$d_uc
  }
)
sigma0_l_values <- lapply(
  X = sigma0_values,
  FUN = function(i) {
    i$l_mat_strict
  }
)
```

### `FitDTVARMxID`

``` r

fit <- FitDTVARMxID(
  data = data,
  observed = c("y1", "y2"),
  id = "id",
  center = TRUE,
  mu_eta_values = mu_eta_values,
  beta_values = beta_values,
  psi_d_values = psi_d_values,
  psi_l_values = psi_l_values,
  theta_d_values = theta_d_values,
  mu0_values = mu0_values,
  sigma0_d_values = sigma0_d_values,
  sigma0_l_values = sigma0_l_values,
  ncores = parallel::detectCores()
)
```

#### Parameter estimates

``` r

head(summary(fit))
#>                             beta_1_1    beta_2_1     beta_1_2  beta_2_2
#> FitDTVARMxID_DTVAR_ID1.Rds 0.7356647 -0.03709508  0.124911294 0.6075660
#> FitDTVARMxID_DTVAR_ID2.Rds 0.3932550  0.33128856  0.070209848 0.6452306
#> FitDTVARMxID_DTVAR_ID3.Rds 0.3667931  0.17275648  0.005867861 0.6412371
#> FitDTVARMxID_DTVAR_ID4.Rds 0.4814453 -0.20665773 -0.083737394 0.1282792
#> FitDTVARMxID_DTVAR_ID5.Rds 0.3894912 -0.09250685  0.029823476 0.5788790
#> FitDTVARMxID_DTVAR_ID6.Rds 0.1466637 -0.04662990  0.044473537 0.3646829
#>                            mu_eta_1_1    mu_eta_2_1   psi_l_2_1   psi_d_1_1
#> FitDTVARMxID_DTVAR_ID1.Rds  0.7429922 -0.5422434386 -0.64137358 -2.09263860
#> FitDTVARMxID_DTVAR_ID2.Rds  2.1089697 -0.9809131058 -0.37478408 -1.74439266
#> FitDTVARMxID_DTVAR_ID3.Rds  0.9790736 -1.7645485016 -0.26088418 -1.21125563
#> FitDTVARMxID_DTVAR_ID4.Rds  1.5471235 -0.7583743612 -0.05025920 -1.63597326
#> FitDTVARMxID_DTVAR_ID5.Rds  0.8140559 -0.0001548687 -0.11296257 -0.87220954
#> FitDTVARMxID_DTVAR_ID6.Rds  0.6719574 -1.2148346783 -0.09648513  0.02601587
#>                              psi_d_2_1 theta_d_1_1 theta_d_2_1
#> FitDTVARMxID_DTVAR_ID1.Rds -2.44661974  -0.3062973  -0.2501242
#> FitDTVARMxID_DTVAR_ID2.Rds -2.01994760  -0.2184480  -0.3630670
#> FitDTVARMxID_DTVAR_ID3.Rds -1.83009313  -0.7312709  -0.5093707
#> FitDTVARMxID_DTVAR_ID4.Rds  0.08626565  -0.3614859 -14.3344963
#> FitDTVARMxID_DTVAR_ID5.Rds -1.77703664  -0.7189124  -0.6065102
#> FitDTVARMxID_DTVAR_ID6.Rds -0.80564615 -17.1757445  -0.8495140
```

#### Proportion of converged cases

``` r

converged(
  fit,
  theta_tol = 0.01,
  prop = TRUE
)
#> [1] 0.904
```

#### Fixed-Effect Meta-Analysis of Measurement Error

When fitting DT-VAR models per person, separating process noise
($`\boldsymbol{\Psi}`$) from measurement error ($`\boldsymbol{\Theta}`$)
can be unstable for some individuals. To stabilize inference, we first
pool the person-level $`\boldsymbol{\Theta}_{i}`$ estimates from only
the converged fits using a fixed-effect meta-analysis. This yields a
high-precision estimate of the common measurement-error covariance that
we will then hold fixed in a second pass of model fitting.

What the code does: - Selects individuals that converged and whose
$`\boldsymbol{\Theta}_i`$ diagonals exceed a small threshold
(`theta_tol`), filtering out near-zero or ill-conditioned solutions. -
Extracts each person’s LDL’ diagonal parameters for
$`\boldsymbol{\Theta}_i`$ and their sampling covariance matrices. -
Computes the inverse-variance-weighted pooled estimate (fixed effect),
returning it on the same LDL’ parameterization used by
[`FitDTVARMxID()`](https://github.com/jeksterslab/fitDTVARMxID/reference/FitDTVARMxID.html).

``` r

library(metaDyn)
fixed_theta <- MetaVARMx(
  fit,
  random = FALSE, # TRUE by default
  effects = FALSE, # TRUE by default
  cov_meas = TRUE, # FALSE by default
  theta_tol = 0.01,
  ncores = parallel::detectCores()
)
```

You can read `summary(fixed_theta)` as providing the pooled (fixed)
measurement-error scale that is common across persons. If individual
instruments truly share the same reliability structure, fixing
$`\boldsymbol{\Theta}`$ to this pooled value improves stability and
often reduces bias in the dynamic parameters.

> **Note:** Fixed-effect pooling assumes a common
> $`\boldsymbol{\Theta}`$ across individuals.

``` r

coef(fixed_theta)
#>  alpha_1_1  alpha_2_1 
#> -0.3799274 -0.3843772
summary(fixed_theta)
#> [1] 0
#> Call:
#> MetaVARMx(object = fit, random = FALSE, effects = FALSE, cov_meas = TRUE, 
#>     theta_tol = 0.01, ncores = parallel::detectCores())
#> 
#> CI type = "normal"
#>                est     se        z p    2.5%   97.5%
#> alpha[1,1] -0.3799 0.0045 -84.0851 0 -0.3888 -0.3711
#> alpha[2,1] -0.3844 0.0043 -88.6964 0 -0.3929 -0.3759
```

``` r

theta_d_values <- coef(fixed_theta)
```

#### Refit the model with fixed measurement error covariance matrix

We refit the individual models using the pooled $`\boldsymbol{\Theta}`$
as a fixed measurement-error covariance matrix.

``` r

fit <- FitDTVARMxID(
  data = data,
  observed = c("y1", "y2"),
  id = "id",
  center = TRUE,
  mu_eta_values = mu_eta_values,
  beta_values = beta_values,
  psi_d_values = psi_d_values,
  psi_l_values = psi_l_values,
  theta_fixed = TRUE,
  theta_d_values = theta_d_values,
  mu0_values = mu0_values,
  sigma0_d_values = sigma0_d_values,
  sigma0_l_values = sigma0_l_values,
  ncores = parallel::detectCores()
)
```

With $`\boldsymbol{\Theta}`$ fixed, the re-estimation focuses on the
dynamic structure ($`\boldsymbol{\mu}`$, $`\boldsymbol{\beta}`$,
$`\boldsymbol{\Psi}`$). In practice, this often increases the proportion
of converged fits and yields more stable cross-lag estimates.

#### Proportion of converged cases

``` r

converged(
  fit,
  prop = TRUE
)
#> Error in `x$output`:
#> ! $ operator is invalid for atomic vectors
```

## Mixed-Effects Meta-Analysis of Person-Specific Dynamics and Means

Having stabilized $`\boldsymbol{\Theta}`$, we synthesize the
person-specific estimates to recover population-level effects,
between-person variability, and systematic covariate-related
differences. We fit a mixed-effects meta-analytic model in which each
individual’s estimate is weighted by its within-person sampling
uncertainty, random effects capture residual heterogeneity across
individuals, and fixed effects quantify the association between
covariate $`X`$ and the person-specific dynamic parameters.

``` r

random <- MetaVARMx(
  fit,
  x = x,
  effects = TRUE,
  set_point = TRUE,
  robust_v = FALSE,
  robust = TRUE,
  ncores = parallel::detectCores()
)
#> Error in `x$output`:
#> ! $ operator is invalid for atomic vectors
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
  person-specific parameters (e.g., $`\mathbb{E}[\boldsymbol{\nu}_i]`$
  and $`\mathbb{E}[\mathrm{vec}(\boldsymbol{\beta}_i)]`$) and fixed
  effects for how $`X`$ shifts these parameters.
- The random part yields between-person covariances
  ($`\boldsymbol{\tau}^2`$) quantifying residual heterogeneity in
  set-point and dynamics beyond what is explained by $`X`$.

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
#> [1]  2.25717566 -2.11275830  0.61362799 -0.00661199 -0.00502123  0.61560881
pop_cov
#>             [,1]        [,2]         [,3]          [,4]          [,5]
#> [1,]  7.87336834  0.88804821  0.273086594  0.0638425889 -0.0915731218
#> [2,]  0.88804821  8.36786378 -0.036344806  0.1331000366 -0.0374675414
#> [3,]  0.27308659 -0.03634481  0.030143625  0.0079252730 -0.0011555838
#> [4,]  0.06384259  0.13310004  0.007925273  0.0139995106 -0.0006635365
#> [5,] -0.09157312 -0.03746754 -0.001155584 -0.0006635365  0.0097715146
#> [6,]  0.09000510 -0.23694437  0.012738301 -0.0012966040  0.0040795653
#>              [,6]
#> [1,]  0.090005095
#> [2,] -0.236944365
#> [3,]  0.012738301
#> [4,] -0.001296604
#> [5,]  0.004079565
#> [6,]  0.027041975
```

## Summary

This vignette demonstrates a two-stage hierarchical estimation approach
for dynamic systems: 1. Individual-level DT-VAR estimation with
stabilized measurement error. 2. Population-level meta-analysis of
person-specific dynamics and means. 3. Estimation and interpretation of
covariate effects, where $`X`$ predicts systematic between-person
differences in dynamics and baseline levels.

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
