# Fit the Discrete-Time Vector Autoregressive Model By ID (Stable Reciprocal Regulation)

## Dynamics Description

The *Stable Reciprocal Regulation* process represents a bivariate
dynamic system in which two latent psychological constructs—such as
positive and negative affect—mutually influence each other over time.
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
component (e.g., positive affect) are naturally counteracted by
adjustments in its counterpart (e.g., negative affect), leading to
emotional homeostasis over time.

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

Let the intercept vector $`\boldsymbol{\alpha}`$ be normally distributed
with the following means
``` math
\begin{equation}
  \left(
    \begin{array}{c}
      0.5 \\
      -0.5 \\
    \end{array}
  \right)
\end{equation}
```
and covariance matrix
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

Let the transition matrix $`\boldsymbol{\beta}`$ be normally distributed
with the following means
``` math
\begin{equation}
  \left(
    \begin{array}{cc}
      0.7 & -0.2 \\
      -0.15 & 0.65 \\
    \end{array}
  \right)
\end{equation}
```
and covariance matrix
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

The `SimAlphaN` and `SimBetaN` functions from the `simStateSpace`
package generate random intercept vectors and transition matrices from
the multivariate normal distribution. Note that the `SimBetaN` function
generates transition matrices that are weakly stationary with an option
to set lower and upper bounds. The person-specific set-point vector
$`\boldsymbol{\mu}_{i}`$ was derived from the generated
$`\boldsymbol{\alpha}_{i}`$ and $`\boldsymbol{\beta}_{i}`$.

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
#> [1]  4.738353 -3.747395
# first sigma0 in the list of length n
sigma0[[1]]
#>            [,1]       [,2]
#> [1,]  1.2855139 -0.9267996
#> [2,] -0.9267996  0.9282909
# first sigma0_l in the list of length n
sigma0_l[[1]] # sigma0_l <- t(chol(sigma0))
#>            [,1]      [,2]
#> [1,]  1.1338050 0.0000000
#> [2,] -0.8174241 0.5100085
# first alpha in the list of length n
alpha[[1]]
#> [1]  0.3319244 -0.1910529
# first beta in the list of length n
beta[[1]]
#>            [,1]       [,2]
#> [1,]  0.6764614 -0.3205201
#> [2,] -0.2392088  0.6465522
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
#> [1]  4.738353 -3.747395
```

### Visualizing the Dynamics Without Process Noise and Measurement Error (n = 5 with Different Initial Condition)

![](fig-vignettes-stable-reciprical-regulation-no-error-1.png)![](fig-vignettes-stable-reciprical-regulation-no-error-2.png)

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
#>   id time       y1        y2
#> 1  1    0 2.196423 -1.922608
#> 2  1    1 1.753266 -3.084149
#> 3  1    2 2.800797 -3.115537
#> 4  1    3 2.654371 -2.559287
#> 5  1    4 3.757502 -1.884908
#> 6  1    5 3.013665 -4.466276
plot(sim, burnin = burnin)
```

![](fig-vignettes-stable-reciprical-regulation-error-1.png)![](fig-vignettes-stable-reciprical-regulation-error-2.png)

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
#>                             beta_1_1    beta_2_1    beta_1_2  beta_2_2
#> FitDTVARMxID_DTVAR_ID1.Rds 0.7368789 -0.16222591 -0.23381631 0.7301704
#> FitDTVARMxID_DTVAR_ID2.Rds 0.7458688  0.08983467 -0.09414717 0.8236796
#> FitDTVARMxID_DTVAR_ID3.Rds 0.8106849 -0.11728965  0.01146235 0.7683006
#> FitDTVARMxID_DTVAR_ID4.Rds 0.7511229 -0.25169156 -0.08074680 0.6060087
#> FitDTVARMxID_DTVAR_ID5.Rds 0.7305870 -0.20750212 -0.13452927 0.6747594
#> FitDTVARMxID_DTVAR_ID6.Rds 0.6585887 -0.23803944 -0.09299216 0.6579883
#>                            mu_eta_1_1 mu_eta_2_1  psi_l_2_1 psi_d_1_1 psi_d_2_1
#> FitDTVARMxID_DTVAR_ID1.Rds   4.756245  -3.779730 -0.3322093 -1.755237 -1.783847
#> FitDTVARMxID_DTVAR_ID2.Rds   4.418347  -2.390498 -0.4988523 -1.698744 -1.926355
#> FitDTVARMxID_DTVAR_ID3.Rds   3.743968  -4.561434 -0.4168239 -1.954509 -1.760115
#> FitDTVARMxID_DTVAR_ID4.Rds   3.756948  -3.135986 -0.2761698 -1.646816 -1.519516
#> FitDTVARMxID_DTVAR_ID5.Rds   2.275892  -1.575123 -0.4945584 -1.478857 -1.586734
#> FitDTVARMxID_DTVAR_ID6.Rds   2.209264  -3.132641 -0.2031644 -1.789331 -1.806059
#>                            theta_d_1_1 theta_d_2_1
#> FitDTVARMxID_DTVAR_ID1.Rds  -0.3486243  -0.5001466
#> FitDTVARMxID_DTVAR_ID2.Rds  -0.2650927  -0.4775165
#> FitDTVARMxID_DTVAR_ID3.Rds  -0.2876378  -0.3585064
#> FitDTVARMxID_DTVAR_ID4.Rds  -0.5770968  -0.5260689
#> FitDTVARMxID_DTVAR_ID5.Rds  -0.5411416  -0.4872107
#> FitDTVARMxID_DTVAR_ID6.Rds  -0.3789378  -0.3863627
```

#### Proportion of converged cases

``` r

converged(
  fit,
  theta_tol = 0.01,
  prop = TRUE
)
#> [1] 0.925
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
#> -0.3907595 -0.3852983
summary(fixed_theta)
#> [1] 0
#> Call:
#> MetaVARMx(object = fit, random = FALSE, effects = FALSE, cov_meas = TRUE, 
#>     theta_tol = 0.01, ncores = parallel::detectCores())
#> 
#> CI type = "normal"
#>                est     se        z p    2.5%   97.5%
#> alpha[1,1] -0.3908 0.0042 -93.4325 0 -0.3990 -0.3826
#> alpha[2,1] -0.3853 0.0043 -89.0446 0 -0.3938 -0.3768
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

## Random-Effects Meta-Analysis of Person-Specific Dynamics and Means

Having stabilized $`\boldsymbol{\Theta}`$, we synthesize the
person-specific estimates to recover population-level effects and their
between-person variability. We use a random-effects model so the pooled
mean reflects both within-person estimation uncertainty and
between-person heterogeneity.

``` r

random <- MetaVARMx(
  fit,
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

- The fixed part of the random-effects model gives pooled means
  $`\boldsymbol{\alpha} = \mathbb{E} \left[ \mathrm{Vec} \left( \boldsymbol{\mu}, \boldsymbol{\beta} \right)  \right]`$.
- The random part yields between-person covariances
  ($`\boldsymbol{\tau}^{2}`$) quantifying heterogeneity in set-point
  ($`\boldsymbol{\mu}`$) and dynamics ($`\boldsymbol{\beta}`$) across
  individuals.

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
#> [1]  4.0081017 -3.1135633  0.6810916 -0.1432002 -0.1942151  0.6367261
pop_cov
#>              [,1]        [,2]          [,3]          [,4]          [,5]
#> [1,] 12.856796905 -5.20249779  0.1964662666  0.0516164907 -0.0883602653
#> [2,] -5.202497788  8.66651926  0.0478426598  0.1419400213 -0.0212600947
#> [3,]  0.196466267  0.04784266  0.0164861144  0.0097360959  0.0002797943
#> [4,]  0.051616491  0.14194002  0.0097360959  0.0144313848 -0.0006153941
#> [5,] -0.088360265 -0.02126009  0.0002797943 -0.0006153941  0.0101332601
#> [6,]  0.003494568 -0.11120175 -0.0001206516  0.0006746266  0.0054812485
#>               [,6]
#> [1,]  0.0034945683
#> [2,] -0.1112017466
#> [3,] -0.0001206516
#> [4,]  0.0006746266
#> [5,]  0.0054812485
#> [6,]  0.0133300763
```

## Summary

This vignette demonstrates a two-stage hierarchical estimation approach
for dynamic systems: 1. individual-level DT-VAR estimation with
stabilized measurement error, and  
2. population-level meta-analysis of person-specific dynamics and means.

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
