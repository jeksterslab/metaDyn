# Meta-Analysis of Continuous-Time Vector Autoregressive Model Estimates

## Dynamics Description

The *Stable Reciprocal Regulation* process is modeled as a
continuous-time bivariate dynamic system in which two latent
psychological constructs (e.g., negative and positive affect) mutually
influence each other through their instantaneous rates of change. The
dynamics are governed by the continuous-time drift matrix whose negative
diagonal elements indicate substantial self-regulation: deviations from
a person’s equilibrium are pulled back toward baseline at a
moderate-to-strong rate. The negative off-diagonal elements reflect
reciprocal inhibition: higher levels of one construct are associated
with an increased instantaneous tendency for the other construct to
decline. Because the cross-effects are modest relative to the
self-regulatory terms, the system exhibits stable, non-oscillatory
relaxation back toward person-specific equilibrium points.

Individuals are allowed to differ in both their self-regulatory rates
and the strength of these antagonistic couplings. Stochastic
disturbances enter through the process-noise covariance, permitting
small (potentially correlated) departures from the deterministic drift,
while measurement errors are assumed to be minimal and comparable across
indicators. Overall, this pattern captures a psychologically plausible
mechanism of reciprocal inhibition in which short-term perturbations in
one component (e.g., negative affect) are naturally counteracted by
compensatory adjustment in the other (e.g., positive affect), supporting
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
  \mathrm{d} \boldsymbol{\eta}_{i, t} = \left( \boldsymbol{\alpha}_{i} + \boldsymbol{\beta}_{i} \boldsymbol{\eta}_{i, t} \right) \mathrm{d} t + \boldsymbol{\Psi}^{\frac{1}{2}} \mathrm{d} \mathbf{W}_{i, t}
\end{equation}
```
where $`\mathrm{d}\boldsymbol{W}_{i, t}`$ is a Wiener process or
Brownian motion, which represents random fluctuations,
$`\boldsymbol{\eta}_{i, t}`$ is a random variable, and
$`\boldsymbol{\beta}_{i}`$, and $`\boldsymbol{\Psi}`$ are model
parameters. Here, $`\boldsymbol{\eta}_{i, t}`$ is a vector of latent
variables at time $`t`$ and individual $`i`$. $`\boldsymbol{\beta}_{i}`$
is a matrix of auto and cross effect coefficients for individual $`i`$,
and $`\boldsymbol{\Psi}`$ the covariance matrix of the process noise.

### Alternative Parameterization

An alternative parameterization of the dynamic structure that directly
estimates the set point vector $`\boldsymbol{\mu}_{i}`$ is given by
``` math
\begin{equation}
  \mathrm{d} \boldsymbol{\eta}_{i, t} = \boldsymbol{\beta}_{i} \left( \boldsymbol{\eta}_{i, t} - \boldsymbol{\mu}_{i} \right) \mathrm{d} t + \boldsymbol{\Psi}^{\frac{1}{2}} \mathrm{d} \mathbf{W}_{i, t} .
\end{equation}
```

Algebraic manipulation of the equation results in the following
``` math
\begin{equation}
  \mathrm{d} \boldsymbol{\eta}_{i, t} =  \left( - \boldsymbol{\beta}_{i} \boldsymbol{\mu}_{i} + \boldsymbol{\beta}_{i} \boldsymbol{\eta}_{i, t} \right) \mathrm{d} t + \boldsymbol{\Psi}^{\frac{1}{2}} \mathrm{d} \mathbf{W}_{i, t}
\end{equation}
```
where we can see that the intercept vector $`\boldsymbol{\alpha}_{i}`$
is implied by $`- \boldsymbol{\beta}_{i} \boldsymbol{\mu}_{i}`$.

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

Let the drift matrix $`\boldsymbol{\beta}`$ be normally distributed with
the following means
``` math
\begin{equation}
  \left(
    \begin{array}{cc}
      -1.2843504 & -0.1792463 \\
      -0.1307004 & -1.3590363 \\
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

The `SimMuN` and `SimPhiN` functions from the `simStateSpace` package
generate random set point vectors and drift matrices from the
multivariate normal distribution. Note that the `SimPhiN` function
generates drift matrices that are weakly stationary with an option to
set lower and upper bounds.

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

\
`n`\
`#> [1] 100`\
`time`\
`#> [1] 1000`\
`# first mu0 in the list of length n`\
`mu0``[[``1``]``]`\
`#> [1] 2.287769 2.606098`\
`# first sigma0 in the list of length n`\
`sigma0``[[``1``]``]`\
`#>           [,1]      [,2]`\
`#> [1,] 0.5037766 0.1478276`\
`#> [2,] 0.1478276 0.5694164`\
`# first sigma0_l in the list of length n`\
`sigma0_l``[[``1``]``]`` ``# sigma0_l <- t(chol(sigma0))`\
`#>           [,1]      [,2]`\
`#> [1,] 0.7097722 0.0000000`\
`#> [2,] 0.2082747 0.7252848`\
`# first alpha in the list of length n`\
`alpha``[[``1``]``]`\
`#> [1] 3.312112 3.842713`\
`# first beta in the list of length n`\
`beta``[[``1``]``]`\
`#>           [,1]      [,2]`\
`#> [1,] -1.235608 -0.186227`\
`#> [2,] -0.169327 -1.325864`\
`# first psi in the list of length n`\
`psi``[[``1``]``]`\
`#>      [,1] [,2]`\
`#> [1,] 1.30 0.57`\
`#> [2,] 0.57 1.56`\
`psi_l``[[``1``]``]`` ``# psi_l <- t(chol(psi))`\
`#>           [,1]     [,2]`\
`#> [1,] 1.1401754 0.000000`\
`#> [2,] 0.4999231 1.144586`\
`# first mu (set point) in the list of length n`\
`mu``[[``1``]``]`\
`#> [1] 2.287769 2.606098`

### Visualizing the Dynamics Without Process Noise and Measurement Error (n = 5 with Different Initial Condition)

![](fig-vignettes-meta-ct-var-no-error-1.png)![](fig-vignettes-meta-ct-var-no-error-2.png)

### Using the `SimSSMLinSDEFixed` Function from the `simStateSpace` Package to Simulate Data

> **Note:** The `SimSSMLinSDEFixed` function uses a different set of
> parameter names. See
> [`help(SimSSMLinSDEFixed)`](https://github.com/jeksterslab/simStateSpace/reference/SimSSMLinSDEFixed.html)
> for more details.

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`simStateSpace`](https://github.com/jeksterslab/simStateSpace)`)`\
`sim`` ``<-`` `[`SimSSMLinSDEIVary`](https://github.com/jeksterslab/simStateSpace/reference/SimSSMLinSDEIVary.html)`(`\
`  n ``=`` ``n``,`\
`  time ``=`` ``time``,`\
`  delta_t ``=`` ``0.1``,`\
`  mu0 ``=`` ``mu0``,`\
`  sigma0_l ``=`` ``sigma0_l``,`\
`  iota ``=`` ``alpha``,`\
`  phi ``=`` ``beta``,`\
`  sigma_l ``=`` ``psi_l``,`\
`  nu ``=`` ``nu``,`\
`  lambda ``=`` ``lambda``,`\
`  theta_l ``=`` ``theta_l`\
`)`\
`data`` ``<-`` `[`as.data.frame`](https://rdrr.io/r/base/as.data.frame.html)`(``sim``)`\
[`head`](https://rdrr.io/r/utils/head.html)`(``data``)`\
`#>   id time       y1       y2`\
`#> 1  1  0.0 2.523315 3.353060`\
`#> 2  1  0.1 2.309397 3.020899`\
`#> 3  1  0.2 2.600318 3.078564`\
`#> 4  1  0.3 2.698363 2.929001`\
`#> 5  1  0.4 2.444623 2.861375`\
`#> 6  1  0.5 3.190602 2.891174`\
[`summary`](https://rdrr.io/r/base/summary.html)`(``data``)`\
`#>        id              time             y1               y2        `\
`#>  Min.   :  1.00   Min.   : 0.00   Min.   :-2.335   Min.   :-1.951  `\
`#>  1st Qu.: 25.75   1st Qu.:24.98   1st Qu.: 2.093   1st Qu.: 1.261  `\
`#>  Median : 50.50   Median :49.95   Median : 2.990   Median : 2.187  `\
`#>  Mean   : 50.50   Mean   :49.95   Mean   : 2.975   Mean   : 2.185  `\
`#>  3rd Qu.: 75.25   3rd Qu.:74.92   3rd Qu.: 3.873   3rd Qu.: 3.081  `\
`#>  Max.   :100.00   Max.   :99.90   Max.   : 8.893   Max.   : 6.970`\
[`plot`](https://rdrr.io/r/graphics/plot.default.html)`(``sim``)`

![](fig-vignettes-meta-ct-var-error-1.png)![](fig-vignettes-meta-ct-var-error-2.png)

## Stage 1: Person-Specific VAR Model

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`OpenMx`](https://openmx.ssri.psu.edu/)`)`\
[`library`](https://rdrr.io/r/base/library.html)`(`[`fitVARMxID`](https://github.com/jeksterslab/fitVARMxID)`)`

The `FitVARMxID` function fits a VAR model on each individual $`i`$.

> **Note:** Consider using the argument `ncores` to use multiple cores
> for parallel processing.

\
`fit`` ``<-`` `[`FitVARMxID`](https://github.com/jeksterslab/fitVARMxID/reference/FitVARMxID.html)`(`\
`  data ``=`` ``data``,`\
`  observed ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"y1"``, ``"y2"``)``,`\
`  id ``=`` ``"id"``,`\
`  ct ``=`` ``TRUE``,`\
`  time ``=`` ``"time"``,`\
`  center ``=`` ``TRUE``,`\
`  ncores ``=`` ``parallel``::`[`detectCores`](https://rdrr.io/r/parallel/detectCores.html)`(``)`\
`)`

\
[`summary`](https://rdrr.io/r/base/summary.html)`(`\
`  ``fit``,`\
`  means ``=`` ``TRUE`\
`)`\
`#> Call:`\
`#> FitVARMxID(data = data, observed = c("y1", "y2"), id = "id", `\
`#>     time = "time", ct = TRUE, center = TRUE, ncores = parallel::detectCores())`\
`#> `\
`#> Convergence:`\
`#> 100.0%`\
`#> `\
`#> Means of the estimated paramaters per individual.`\
`#>   mu[1,1]   mu[2,1] beta[1,1] beta[2,1] beta[1,2] beta[2,2]  psi[1,1]  psi[2,1] `\
`#>    2.9743    2.1847   -1.3528   -0.1673   -0.1959   -1.4373    1.3040    0.5685 `\
`#>  psi[2,2] `\
`#>    1.5566`

## Stage 2: Random-Effects Meta-Analysis of Person-Specific Set Points and Dynamics

We synthesize the person-specific estimates to recover population-level
effects and their between-person variability. We use a random-effects
model so the pooled mean reflects both within-person estimation
uncertainty and between-person heterogeneity.

All available parameters are meta-analyzed by default. Setting
`cov_dyn = FALSE`, meta-analyzes only the set points and drift matrix.
Setting `tau_sqr_l_free`, such that covariances between `mu` and `beta`
are constained to zero, simplifies the random effects.

\
`tau_sqr_l_free`` ``<-`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(`\
`  data ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(`\
`    ``FALSE``, ``FALSE``, ``FALSE``, ``FALSE``, ``FALSE``, ``FALSE``,`\
`    ``TRUE``, ``FALSE``, ``FALSE``, ``FALSE``, ``FALSE``, ``FALSE``,`\
`    ``FALSE``, ``FALSE``, ``FALSE``, ``FALSE``, ``FALSE``, ``FALSE``,`\
`    ``FALSE``, ``FALSE``, ``TRUE``, ``FALSE``, ``FALSE``, ``FALSE``,`\
`    ``FALSE``, ``FALSE``, ``TRUE``, ``TRUE``, ``FALSE``, ``FALSE``,`\
`    ``FALSE``, ``FALSE``, ``TRUE``, ``TRUE``, ``TRUE``, ``FALSE`\
`  ``)``,`\
`  byrow ``=`` ``TRUE``,`\
`  nrow ``=`` ``6``,`\
`  ncol ``=`` ``6`\
`)`\
`# FALSE values in tau_sqr_l_free will be treated as fixed values.`\
`# TRUE values in tau_sqr_l_free will be treated as starting values.`\
`tau_sqr_l_values`` ``<-`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(`\
`  data ``=`` ``0``,`\
`  nrow ``=`` ``6``,`\
`  ncol ``=`` ``6`\
`)`

> **Note:** Consider using the argument `ncores` to use multiple cores
> for parallel processing.

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`metaDyn`](https://github.com/jeksterslab/metaDyn)`)`\
`metavar`` ``<-`` `[`MetaVARMx`](https://github.com/jeksterslab/metaDyn/reference/MetaVARMx.md)`(`\
`  object ``=`` ``fit``,`\
`  tau_sqr_l_free ``=`` ``tau_sqr_l_free``,`\
`  tau_sqr_l_values ``=`` ``tau_sqr_l_values``,`\
`  cov_dyn ``=`` ``FALSE``,`\
`  ncores ``=`` ``parallel``::`[`detectCores`](https://rdrr.io/r/parallel/detectCores.html)`(``)`\
`)`

\
[`summary`](https://rdrr.io/r/base/summary.html)`(``metavar``)`\
`#> Call:`\
`#> MetaVARMx(object = fit, tau_sqr_l_free = tau_sqr_l_free, tau_sqr_l_values = tau_sqr_l_values, `\
`#>     cov_dyn = FALSE, ncores = parallel::detectCores())`\
`#> `\
`#> Status code:`\
`#> 0`\
`#> `\
`#> Confidence intervals type:`\
`#> Wald`\
`#> `\
`#> Sampling covariance matrix type:`\
`#> Normal`\
`#> `\
`#>                  est     se         z      p    2.5%   97.5%`\
`#> alpha[1,1]    2.9742 0.1120   26.5554 0.0000  2.7547  3.1937`\
`#> alpha[2,1]    2.1845 0.1062   20.5782 0.0000  1.9765  2.3926`\
`#> alpha[3,1]   -1.3129 0.0238  -55.2177 0.0000 -1.3595 -1.2663`\
`#> alpha[4,1]   -0.1670 0.0225   -7.4214 0.0000 -0.2111 -0.1229`\
`#> alpha[5,1]   -0.1943 0.0186  -10.4605 0.0000 -0.2307 -0.1579`\
`#> alpha[6,1]   -1.3984 0.0251  -55.8097 0.0000 -1.4475 -1.3493`\
`#> tau_sqr[1,1]  1.2471 0.1775    7.0249 0.0000  0.8991  1.5950`\
`#> tau_sqr[2,1]  0.5434 0.1307    4.1583 0.0000  0.2873  0.7996`\
`#> tau_sqr[2,2]  1.1190 0.1593    7.0230 0.0000  0.8067  1.4313`\
`#> tau_sqr[3,3]  0.0209 0.0072    2.8955 0.0038  0.0068  0.0351`\
`#> tau_sqr[4,3]  0.0013 0.0048    0.2591 0.7956 -0.0082  0.0108`\
`#> tau_sqr[5,3]  0.0005 0.0041    0.1156 0.9080 -0.0075  0.0084`\
`#> tau_sqr[6,3] -0.0046 0.0058   -0.7896 0.4297 -0.0159  0.0068`\
`#> tau_sqr[4,4]  0.0088 0.0063    1.3975 0.1623 -0.0036  0.0212`\
`#> tau_sqr[5,4]  0.0045 0.0032    1.4277 0.1534 -0.0017  0.0107`\
`#> tau_sqr[6,4] -0.0038 0.0056   -0.6870 0.4921 -0.0147  0.0071`\
`#> tau_sqr[5,5]  0.0035 0.0036    0.9686 0.3327 -0.0036  0.0106`\
`#> tau_sqr[6,5]  0.0032 0.0042    0.7704 0.4411 -0.0050  0.0114`\
`#> tau_sqr[6,6]  0.0251 0.0079    3.1573 0.0016  0.0095  0.0406`\
`#> i_sqr[1,1]    0.9948 0.0007 1339.2107 0.0000  0.9933  0.9962`\
`#> i_sqr[2,1]    0.9937 0.0009 1110.1538 0.0000  0.9919  0.9954`\
`#> i_sqr[3,1]    0.3831 0.0816    4.6937 0.0000  0.2231  0.5430`\
`#> i_sqr[4,1]    0.1798 0.1055    1.7040 0.0884 -0.0270  0.3867`\
`#> i_sqr[5,1]    0.1045 0.0966    1.0816 0.2794 -0.0848  0.2938`\
`#> i_sqr[6,1]    0.4128 0.0768    5.3777 0.0000  0.2623  0.5633`

### Normal Theory Confidence Intervals

\
[`confint`](https://rdrr.io/r/stats/confint.html)`(``metavar``, level ``=`` ``0.95``)`\
`#>                     2.5 %       97.5 %`\
`#> alpha[1,1]    2.754699455  3.193733151`\
`#> alpha[2,1]    1.976476141  2.392609090`\
`#> alpha[3,1]   -1.359465091 -1.266264257`\
`#> alpha[4,1]   -0.211107619 -0.122898296`\
`#> alpha[5,1]   -0.230722568 -0.157905851`\
`#> alpha[6,1]   -1.447482343 -1.349264193`\
`#> tau_sqr[1,1]  0.899138573  1.595018070`\
`#> tau_sqr[2,1]  0.287298419  0.799589635`\
`#> tau_sqr[2,2]  0.806734694  1.431327291`\
`#> tau_sqr[3,3]  0.006754915  0.035058325`\
`#> tau_sqr[4,3] -0.008244302  0.010755534`\
`#> tau_sqr[5,3] -0.007500712  0.008440829`\
`#> tau_sqr[6,3] -0.015888896  0.006762876`\
`#> tau_sqr[4,4] -0.003552467  0.021204845`\
`#> tau_sqr[5,4] -0.001688715  0.010749140`\
`#> tau_sqr[6,4] -0.014716571  0.007077148`\
`#> tau_sqr[5,5] -0.003581436  0.010580102`\
`#> tau_sqr[6,5] -0.004981648  0.011433757`\
`#> tau_sqr[6,6]  0.009508950  0.040640207`\
`#> i_sqr[1,1]    0.993298532  0.996210221`\
`#> i_sqr[2,1]    0.991919382  0.995428022`\
`#> i_sqr[3,1]    0.223117850  0.543046403`\
`#> i_sqr[4,1]   -0.027017517  0.386716116`\
`#> i_sqr[5,1]   -0.084848123  0.293830188`\
`#> i_sqr[6,1]    0.262349390  0.563251337`

\
[`confint`](https://rdrr.io/r/stats/confint.html)`(``metavar``, level ``=`` ``0.99``)`\
`#>                     0.5 %      99.5 %`\
`#> alpha[1,1]    2.685722262  3.26271034`\
`#> alpha[2,1]    1.911096916  2.45798832`\
`#> alpha[3,1]   -1.374108003 -1.25162135`\
`#> alpha[4,1]   -0.224966308 -0.10903961`\
`#> alpha[5,1]   -0.242162904 -0.14646552`\
`#> alpha[6,1]   -1.462913533 -1.33383300`\
`#> tau_sqr[1,1]  0.789807980  1.70434866`\
`#> tau_sqr[2,1]  0.206811636  0.88007642`\
`#> tau_sqr[2,2]  0.708604085  1.52945790`\
`#> tau_sqr[3,3]  0.002308127  0.03950511`\
`#> tau_sqr[4,3] -0.011229393  0.01374063`\
`#> tau_sqr[5,3] -0.010005310  0.01094543`\
`#> tau_sqr[6,3] -0.019447748  0.01032173`\
`#> tau_sqr[4,4] -0.007442123  0.02509450`\
`#> tau_sqr[5,4] -0.003642843  0.01270327`\
`#> tau_sqr[6,4] -0.018140613  0.01050119`\
`#> tau_sqr[5,5] -0.005806375  0.01280504`\
`#> tau_sqr[6,5] -0.007560696  0.01401280`\
`#> tau_sqr[6,6]  0.004617875  0.04553128`\
`#> i_sqr[1,1]    0.992841072  0.99666768`\
`#> i_sqr[2,1]    0.991368135  0.99597927`\
`#> i_sqr[3,1]    0.172853432  0.59331082`\
`#> i_sqr[4,1]   -0.092019782  0.45171838`\
`#> i_sqr[5,1]   -0.144342798  0.35332486`\
`#> i_sqr[6,1]    0.215074268  0.61052646`

### Robust Confidence Intervals

\
[`confint`](https://rdrr.io/r/stats/confint.html)`(``metavar``, level ``=`` ``0.95``, robust ``=`` ``TRUE``)`\
`#>                     2.5 %       97.5 %`\
`#> alpha[1,1]    2.753596981  3.194835624`\
`#> alpha[2,1]    1.975423783  2.393661448`\
`#> alpha[3,1]   -1.358269316 -1.267460032`\
`#> alpha[4,1]   -0.211821215 -0.122184700`\
`#> alpha[5,1]   -0.230867134 -0.157761285`\
`#> alpha[6,1]   -1.445594684 -1.351151852`\
`#> tau_sqr[1,1]  0.908547294  1.585609348`\
`#> tau_sqr[2,1]  0.284671621  0.802216434`\
`#> tau_sqr[2,2]  0.848745083  1.389316902`\
`#> tau_sqr[3,3]  0.006112196  0.035701045`\
`#> tau_sqr[4,3] -0.007479194  0.009990427`\
`#> tau_sqr[5,3] -0.007526798  0.008466915`\
`#> tau_sqr[6,3] -0.015336194  0.006210173`\
`#> tau_sqr[4,4] -0.002288021  0.019940398`\
`#> tau_sqr[5,4] -0.001041266  0.010101691`\
`#> tau_sqr[6,4] -0.012814010  0.005174588`\
`#> tau_sqr[5,5] -0.002924336  0.009923001`\
`#> tau_sqr[6,5] -0.004272275  0.010724383`\
`#> tau_sqr[6,6]  0.005247552  0.044901605`\
`#> i_sqr[1,1]    0.993337900  0.996170853`\
`#> i_sqr[2,1]    0.992155375  0.995192029`\
`#> i_sqr[3,1]    0.215852848  0.550311405`\
`#> i_sqr[4,1]   -0.005886582  0.365585181`\
`#> i_sqr[5,1]   -0.067277306  0.276259371`\
`#> i_sqr[6,1]    0.221136097  0.604464630`

\
[`confint`](https://rdrr.io/r/stats/confint.html)`(``metavar``, level ``=`` ``0.99``, robust ``=`` ``TRUE``)`\
`#>                      0.5 %       99.5 %`\
`#> alpha[1,1]    2.6842733665  3.264159239`\
`#> alpha[2,1]    1.9097138830  2.459371348`\
`#> alpha[3,1]   -1.3725364890 -1.253192859`\
`#> alpha[4,1]   -0.2259041322 -0.108101783`\
`#> alpha[5,1]   -0.2423528958 -0.146275524`\
`#> alpha[6,1]   -1.4604327287 -1.336313807`\
`#> tau_sqr[1,1]  0.8021731363  1.691983506`\
`#> tau_sqr[2,1]  0.2033594379  0.883528617`\
`#> tau_sqr[2,2]  0.7638150951  1.474246890`\
`#> tau_sqr[3,3]  0.0014634508  0.040349790`\
`#> tau_sqr[4,3] -0.0102238708  0.012735103`\
`#> tau_sqr[5,3] -0.0100395919  0.010979709`\
`#> tau_sqr[6,3] -0.0187213731  0.009595352`\
`#> tau_sqr[4,4] -0.0057803589  0.023432736`\
`#> tau_sqr[5,4] -0.0027919509  0.011852376`\
`#> tau_sqr[6,4] -0.0156402239  0.008000801`\
`#> tau_sqr[5,5] -0.0049427988  0.011941464`\
`#> tau_sqr[6,5] -0.0066284206  0.013080529`\
`#> tau_sqr[6,6] -0.0009825516  0.051131708`\
`#> i_sqr[1,1]    0.9928928104  0.996615943`\
`#> i_sqr[2,1]    0.9916782825  0.995669122`\
`#> i_sqr[3,1]    0.1633056003  0.602858653`\
`#> i_sqr[4,1]   -0.0642490257  0.423947625`\
`#> i_sqr[5,1]   -0.1212508298  0.330232894`\
`#> i_sqr[6,1]    0.1609108194  0.664689907`

- The fixed part of the random-effects model gives pooled means
  $`\boldsymbol{\alpha} = \mathbb{E} \left[ \mathrm{Vec} \left( \boldsymbol{\mu}, \boldsymbol{\beta} \right)  \right]`$.
- The random part yields between-person covariances
  ($`\boldsymbol{\tau}^{2}`$) quantifying heterogeneity in set point
  ($`\boldsymbol{\mu}`$) and dynamics ($`\boldsymbol{\beta}`$) across
  individuals.

\
`means`` ``<-`` `[`extract`](https://github.com/jeksterslab/metaDyn/reference/extract.md)`(``metavar``, what ``=`` ``"alpha"``)`\
`means`\
`#>         alpha`\
`#> y1  2.9742163`\
`#> y2  2.1845426`\
`#> y3 -1.3128647`\
`#> y4 -0.1670030`\
`#> y5 -0.1943142`\
`#> y6 -1.3983733`\
`covariances`` ``<-`` `[`extract`](https://github.com/jeksterslab/metaDyn/reference/extract.md)`(``metavar``, what ``=`` ``"tau_sqr"``)`\
`covariances`\
`#>          y1       y2            y3           y4           y5           y6`\
`#> y1 1.247078 0.543444  0.0000000000  0.000000000 0.0000000000  0.000000000`\
`#> y2 0.543444 1.119031  0.0000000000  0.000000000 0.0000000000  0.000000000`\
`#> y3 0.000000 0.000000  0.0209066202  0.001255616 0.0004700585 -0.004563010`\
`#> y4 0.000000 0.000000  0.0012556161  0.008826189 0.0045302125 -0.003819711`\
`#> y5 0.000000 0.000000  0.0004700585  0.004530213 0.0034993328  0.003226054`\
`#> y6 0.000000 0.000000 -0.0045630104 -0.003819711 0.0032260541  0.025074578`

## Summary

This vignette demonstrates a two-stage hierarchical estimation approach
for dynamic systems: 1. individual-level CT-VAR estimation, and\
2. population-level meta-analysis of person-specific set points and
dynamics.

## References

Cheung, M. W.-L. (2015). *Meta-analysis: A structural equation modeling
approach*. Wiley. <https://doi.org/10.1002/9781118957813>

R Core Team. (2026). *R: A language and environment for statistical
computing*. R Foundation for Statistical Computing.
<https://www.R-project.org/>
