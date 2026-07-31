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

`n`` ``#> [1] 100`` ``time`` ``#> [1] 1000`` ``# first mu0 in the list of length n`` ``mu0``[[``1``]``]`` ``#> [1] 2.287769 2.606098`` ``# first sigma0 in the list of length n`` ``sigma0``[[``1``]``]`` ``#> [,1] [,2]`` ``#> [1,] 1.4403549 0.5646226`` ``#> [2,] 0.5646226 1.6885649`` ``# first sigma0_l in the list of length n`` ``sigma0_l``[[``1``]``]`` ``# sigma0_l <- t(chol(sigma0))`` ``#> [,1] [,2]`` ``#> [1,] 1.2001479 0.000000`` ``#> [2,] 0.4704608 1.211293`` ``# first alpha in the list of length n`` ``alpha``[[``1``]``]`` ``#> [1] 1.678968 2.010503`` ``# first beta in the list of length n`` ``beta``[[``1``]``]`` ``#> [,1] [,2]`` ``#> [1,] 0.32874206 -0.05498069`` ``#> [2,] -0.07362662 0.29317209`` ``# first psi in the list of length n`` ``psi``[[``1``]``]`` ``#> [,1] [,2]`` ``#> [1,] 1.30 0.57`` ``#> [2,] 0.57 1.56`` ``psi_l``[[``1``]``]`` ``# psi_l <- t(chol(psi))`` ``#> [,1] [,2]`` ``#> [1,] 1.1401754 0.000000`` ``#> [2,] 0.4999231 1.144586`` ``# first mu (set point) in the list of length n`` ``mu``[[``1``]``]`` ``#> [1] 2.287769 2.606098`

### Visualizing the Dynamics Without Process Noise and Measurement Error (n = 5 with Different Initial Condition)

![](fig-vignettes-meta-dt-var-no-error-1.png)![](fig-vignettes-meta-dt-var-no-error-2.png)

### Using the `SimSSMVARIVary` Function from the `simStateSpace` Package to Simulate Data

[`library`](https://rdrr.io/r/base/library.html)`(`[`simStateSpace`](https://github.com/jeksterslab/simStateSpace)`)`` ``sim`` ``<-`` `[`SimSSMVARIVary`](https://github.com/jeksterslab/simStateSpace/reference/SimSSMVARIVary.html)`(`` `` n ``=`` ``n``,`` `` time ``=`` ``time``,`` `` mu0 ``=`` ``mu0``,`` `` sigma0_l ``=`` ``sigma0_l``,`` `` alpha ``=`` ``alpha``,`` `` beta ``=`` ``beta``,`` `` psi_l ``=`` ``psi_l`` ``)`` ``data`` ``<-`` `[`as.data.frame`](https://rdrr.io/r/base/as.data.frame.html)`(``sim``)`` `[`head`](https://rdrr.io/r/utils/head.html)`(``data``)`` ``#> id time y1 y2`` ``#> 1 1 0 2.2494783 1.541912`` ``#> 2 1 1 2.8646254 2.346091`` ``#> 3 1 2 0.8218736 3.348071`` ``#> 4 1 3 2.4900816 2.282659`` ``#> 5 1 4 2.7074414 3.944741`` ``#> 6 1 5 2.9755157 2.667683`` `[`plot`](https://rdrr.io/r/graphics/plot.default.html)`(``sim``)`

![](fig-vignettes-meta-dt-var-error-1.png)![](fig-vignettes-meta-dt-var-error-2.png)

## Stage 1: Person-Specific VAR Model

[`library`](https://rdrr.io/r/base/library.html)`(`[`OpenMx`](https://openmx.ssri.psu.edu/)`)`` `[`library`](https://rdrr.io/r/base/library.html)`(`[`fitVARMxID`](https://github.com/jeksterslab/fitVARMxID)`)`

The `FitVARMxID` function fits a VAR model on each individual $`i`$.

> **Note:** Consider using the argument `ncores` to use multiple cores
> for parallel processing.

`fit`` ``<-`` `[`FitVARMxID`](https://github.com/jeksterslab/fitVARMxID/reference/FitVARMxID.html)`(`` `` data ``=`` ``data``,`` `` observed ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"y1"``, ``"y2"``)``,`` `` id ``=`` ``"id"``,`` `` center ``=`` ``TRUE``,`` `` ncores ``=`` ``parallel``::`[`detectCores`](https://rdrr.io/r/parallel/detectCores.html)`(``)`` ``)`

[`summary`](https://rdrr.io/r/base/summary.html)`(`` `` ``fit``,`` `` means ``=`` ``TRUE`` ``)`` ``#> Call:`` ``#> FitVARMxID(data = data, observed = c("y1", "y2"), id = "id", `` ``#> center = TRUE, ncores = parallel::detectCores())`` ``#> `` ``#> Convergence:`` ``#> 100.0%`` ``#> `` ``#> Means of the estimated paramaters per individual.`` ``#> mu[1,1] mu[2,1] beta[1,1] beta[2,1] beta[1,2] beta[2,2] psi[1,1] psi[2,1] `` ``#> 2.9733 2.1842 0.2711 -0.0659 -0.0543 0.2394 1.2932 0.5612 `` ``#> psi[2,2] `` ``#> 1.5390`

## Stage 2: Random-Effects Meta-Analysis of Person-Specific Set Points and Dynamics

We synthesize the person-specific estimates to recover population-level
effects and their between-person variability. We use a random-effects
model so the pooled mean reflects both within-person estimation
uncertainty and between-person heterogeneity.

All available parameters are meta-analyzed by default. Setting
`cov_dyn = FALSE`, meta-analyzes only the set points and transition
matrix. Setting `tau_sqr_l_free`, such that covariances between `mu` and
`beta` are constained to zero, simplifies the random effects.

`tau_sqr_l_free`` ``<-`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(`` `` data ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(`` `` ``FALSE``, ``FALSE``, ``FALSE``, ``FALSE``, ``FALSE``, ``FALSE``,`` `` ``TRUE``, ``FALSE``, ``FALSE``, ``FALSE``, ``FALSE``, ``FALSE``,`` `` ``FALSE``, ``FALSE``, ``FALSE``, ``FALSE``, ``FALSE``, ``FALSE``,`` `` ``FALSE``, ``FALSE``, ``TRUE``, ``FALSE``, ``FALSE``, ``FALSE``,`` `` ``FALSE``, ``FALSE``, ``TRUE``, ``TRUE``, ``FALSE``, ``FALSE``,`` `` ``FALSE``, ``FALSE``, ``TRUE``, ``TRUE``, ``TRUE``, ``FALSE`` `` ``)``,`` `` byrow ``=`` ``TRUE``,`` `` nrow ``=`` ``6``,`` `` ncol ``=`` ``6`` ``)`` ``# FALSE values in tau_sqr_l_free will be treated as fixed values.`` ``# TRUE values in tau_sqr_l_free will be treated as starting values.`` ``tau_sqr_l_values`` ``<-`` `[`matrix`](https://rdrr.io/r/base/matrix.html)`(`` `` data ``=`` ``0``,`` `` nrow ``=`` ``6``,`` `` ncol ``=`` ``6`` ``)`

> **Note:** Consider using the argument `ncores` to use multiple cores
> for parallel processing.

[`library`](https://rdrr.io/r/base/library.html)`(`[`metaDyn`](https://github.com/jeksterslab/metaDyn)`)`` ``metavar`` ``<-`` `[`MetaVARMx`](https://github.com/jeksterslab/metaDyn/reference/MetaVARMx.md)`(`` `` object ``=`` ``fit``,`` `` tau_sqr_l_free ``=`` ``tau_sqr_l_free``,`` `` tau_sqr_l_values ``=`` ``tau_sqr_l_values``,`` `` cov_dyn ``=`` ``FALSE``,`` `` ncores ``=`` ``parallel``::`[`detectCores`](https://rdrr.io/r/parallel/detectCores.html)`(``)`` ``)`

[`summary`](https://rdrr.io/r/base/summary.html)`(``metavar``)`` ``#> Call:`` ``#> MetaVARMx(object = fit, tau_sqr_l_free = tau_sqr_l_free, tau_sqr_l_values = tau_sqr_l_values, `` ``#> cov_dyn = FALSE, ncores = parallel::detectCores())`` ``#> `` ``#> Status code:`` ``#> 0`` ``#> `` ``#> Confidence intervals type:`` ``#> Wald`` ``#> `` ``#> Sampling covariance matrix type:`` ``#> Normal`` ``#> `` ``#> est se z p 2.5% 97.5%`` ``#> alpha[1,1] 2.9733 0.1115 26.6781 0.0000 2.7549 3.1918`` ``#> alpha[2,1] 2.1844 0.1052 20.7594 0.0000 1.9781 2.3906`` ``#> alpha[3,1] 0.2720 0.0122 22.3157 0.0000 0.2481 0.2959`` ``#> alpha[4,1] -0.0659 0.0106 -6.2428 0.0000 -0.0866 -0.0452`` ``#> alpha[5,1] -0.0545 0.0047 -11.6738 0.0000 -0.0636 -0.0453`` ``#> alpha[6,1] 0.2401 0.0150 16.0294 0.0000 0.2107 0.2694`` ``#> tau_sqr[1,1] 1.2396 0.1757 7.0558 0.0000 0.8952 1.5839`` ``#> tau_sqr[2,1] 0.5389 0.1291 4.1756 0.0000 0.2859 0.7918`` ``#> tau_sqr[2,2] 1.1040 0.1566 7.0513 0.0000 0.7971 1.4108`` ``#> tau_sqr[3,3] 0.0138 0.0021 6.6038 0.0000 0.0097 0.0179`` ``#> tau_sqr[4,3] 0.0021 0.0013 1.5812 0.1138 -0.0005 0.0046`` ``#> tau_sqr[5,3] 0.0013 0.0006 2.3450 0.0190 0.0002 0.0025`` ``#> tau_sqr[6,3] 0.0037 0.0019 2.0081 0.0446 0.0001 0.0074`` ``#> tau_sqr[4,4] 0.0099 0.0016 6.2655 0.0000 0.0068 0.0130`` ``#> tau_sqr[5,4] 0.0007 0.0005 1.4928 0.1355 -0.0002 0.0017`` ``#> tau_sqr[6,4] 0.0053 0.0017 3.2207 0.0013 0.0021 0.0086`` ``#> tau_sqr[5,5] 0.0013 0.0003 4.2040 0.0000 0.0007 0.0019`` ``#> tau_sqr[6,5] 0.0023 0.0007 3.1303 0.0017 0.0008 0.0037`` ``#> tau_sqr[6,6] 0.0214 0.0032 6.7778 0.0000 0.0152 0.0276`` ``#> i_sqr[1,1] 0.9982 0.0003 3860.4659 0.0000 0.9977 0.9987`` ``#> i_sqr[2,1] 0.9978 0.0003 3176.9731 0.0000 0.9972 0.9984`` ``#> i_sqr[3,1] 0.9293 0.0100 93.3767 0.0000 0.9098 0.9488`` ``#> i_sqr[4,1] 0.8879 0.0159 55.9006 0.0000 0.8568 0.9190`` ``#> i_sqr[5,1] 0.5975 0.0572 10.4450 0.0000 0.4854 0.7096`` ``#> i_sqr[6,1] 0.9537 0.0065 146.4235 0.0000 0.9409 0.9665`

### Normal Theory Confidence Intervals

[`confint`](https://rdrr.io/r/stats/confint.html)`(``metavar``, level ``=`` ``0.95``)`` ``#> 2.5 % 97.5 %`` ``#> alpha[1,1] 2.754884e+00 3.191767759`` ``#> alpha[2,1] 1.978127e+00 2.390593407`` ``#> alpha[3,1] 2.481421e-01 0.295927240`` ``#> alpha[4,1] -8.662881e-02 -0.045230631`` ``#> alpha[5,1] -6.361186e-02 -0.045322469`` ``#> alpha[6,1] 2.107037e-01 0.269408542`` ``#> tau_sqr[1,1] 8.952314e-01 1.583880798`` ``#> tau_sqr[2,1] 2.859450e-01 0.791843210`` ``#> tau_sqr[2,2] 7.971012e-01 1.410801329`` ``#> tau_sqr[3,3] 9.704781e-03 0.017896798`` ``#> tau_sqr[4,3] -4.944664e-04 0.004622909`` ``#> tau_sqr[5,3] 2.204776e-04 0.002465307`` ``#> tau_sqr[6,3] 8.937825e-05 0.007368280`` ``#> tau_sqr[4,4] 6.796436e-03 0.012984107`` ``#> tau_sqr[5,4] -2.307973e-04 0.001705878`` ``#> tau_sqr[6,4] 2.085967e-03 0.008571952`` ``#> tau_sqr[5,5] 6.881521e-04 0.001890228`` ``#> tau_sqr[6,5] 8.490762e-04 0.003692964`` ``#> tau_sqr[6,6] 1.519139e-02 0.027551412`` ``#> i_sqr[1,1] 9.976655e-01 0.998679032`` ``#> i_sqr[2,1] 9.971649e-01 0.998395994`` ``#> i_sqr[3,1] 9.097638e-01 0.948774233`` ``#> i_sqr[4,1] 8.567711e-01 0.919033616`` ``#> i_sqr[5,1] 4.853538e-01 0.709578318`` ``#> i_sqr[6,1] 9.409421e-01 0.966474018`

[`confint`](https://rdrr.io/r/stats/confint.html)`(``metavar``, level ``=`` ``0.99``)`` ``#> 0.5 % 99.5 %`` ``#> alpha[1,1] 2.6862443912 3.260407199`` ``#> alpha[2,1] 1.9133242933 2.455396513`` ``#> alpha[3,1] 0.2406345579 0.303434823`` ``#> alpha[4,1] -0.0931329360 -0.038726505`` ``#> alpha[5,1] -0.0664853281 -0.042448997`` ``#> alpha[6,1] 0.2014804790 0.278631744`` ``#> tau_sqr[1,1] 0.7870367406 1.692075459`` ``#> tau_sqr[2,1] 0.2064626183 0.871325581`` ``#> tau_sqr[2,2] 0.7006819564 1.507220602`` ``#> tau_sqr[3,3] 0.0084177216 0.019183858`` ``#> tau_sqr[4,3] -0.0012984644 0.005426907`` ``#> tau_sqr[5,3] -0.0001322108 0.002817996`` ``#> tau_sqr[6,3] -0.0010542202 0.008511879`` ``#> tau_sqr[4,4] 0.0058242829 0.013956261`` ``#> tau_sqr[5,4] -0.0005350711 0.002010152`` ``#> tau_sqr[6,4] 0.0010669444 0.009590974`` ``#> tau_sqr[5,5] 0.0004992923 0.002079088`` ``#> tau_sqr[6,5] 0.0004022691 0.004139771`` ``#> tau_sqr[6,6] 0.0132494888 0.029493313`` ``#> i_sqr[1,1] 0.9975062457 0.998838272`` ``#> i_sqr[2,1] 0.9969714541 0.998589417`` ``#> i_sqr[3,1] 0.9036347907 0.954903220`` ``#> i_sqr[4,1] 0.8469889749 0.928815763`` ``#> i_sqr[5,1] 0.4501256134 0.744806538`` ``#> i_sqr[6,1] 0.9369307974 0.970485366`

### Robust Confidence Intervals

[`confint`](https://rdrr.io/r/stats/confint.html)`(``metavar``, level ``=`` ``0.95``, robust ``=`` ``TRUE``)`` ``#> 2.5 % 97.5 %`` ``#> alpha[1,1] 2.7537802400 3.192871350`` ``#> alpha[2,1] 1.9770925393 2.391628267`` ``#> alpha[3,1] 0.2480753691 0.295994012`` ``#> alpha[4,1] -0.0867440516 -0.045115390`` ``#> alpha[5,1] -0.0636813652 -0.045252960`` ``#> alpha[6,1] 0.2105511451 0.269561078`` ``#> tau_sqr[1,1] 0.9013921845 1.577720015`` ``#> tau_sqr[2,1] 0.2794401393 0.798348060`` ``#> tau_sqr[2,2] 0.8278197995 1.380082759`` ``#> tau_sqr[3,3] 0.0091150765 0.018486503`` ``#> tau_sqr[4,3] -0.0002357698 0.004364212`` ``#> tau_sqr[5,3] 0.0003160367 0.002369748`` ``#> tau_sqr[6,3] -0.0005935769 0.008051236`` ``#> tau_sqr[4,4] 0.0072590988 0.012521445`` ``#> tau_sqr[5,4] -0.0001960241 0.001671105`` ``#> tau_sqr[6,4] 0.0023412891 0.008316630`` ``#> tau_sqr[5,5] 0.0007204161 0.001857964`` ``#> tau_sqr[6,5] 0.0011175122 0.003424528`` ``#> tau_sqr[6,6] 0.0134984520 0.029244349`` ``#> i_sqr[1,1] 0.9976745528 0.998669965`` ``#> i_sqr[2,1] 0.9972265000 0.998334371`` ``#> i_sqr[3,1] 0.9069556006 0.951582410`` ``#> i_sqr[4,1] 0.8614265969 0.914378141`` ``#> i_sqr[5,1] 0.4913722361 0.703559915`` ``#> i_sqr[6,1] 0.9374451947 0.969970968`

[`confint`](https://rdrr.io/r/stats/confint.html)`(``metavar``, level ``=`` ``0.99``, robust ``=`` ``TRUE``)`` ``#> 0.5 % 99.5 %`` ``#> alpha[1,1] 2.684794e+00 3.261857563`` ``#> alpha[2,1] 1.911964e+00 2.456756550`` ``#> alpha[3,1] 2.405468e-01 0.303522577`` ``#> alpha[4,1] -9.328439e-02 -0.038575053`` ``#> alpha[5,1] -6.657668e-02 -0.042357648`` ``#> alpha[6,1] 2.012800e-01 0.278832210`` ``#> tau_sqr[1,1] 7.951334e-01 1.683978819`` ``#> tau_sqr[2,1] 1.979138e-01 0.879874402`` ``#> tau_sqr[2,2] 7.410530e-01 1.466849558`` ``#> tau_sqr[3,3] 7.642719e-03 0.019958861`` ``#> tau_sqr[4,3] -9.584793e-04 0.005086922`` ``#> tau_sqr[5,3] -6.624804e-06 0.002692410`` ``#> tau_sqr[6,3] -1.951775e-03 0.009409434`` ``#> tau_sqr[4,4] 6.432324e-03 0.013348219`` ``#> tau_sqr[5,4] -4.893714e-04 0.001964452`` ``#> tau_sqr[6,4] 1.402495e-03 0.009255424`` ``#> tau_sqr[5,5] 5.416943e-04 0.002036686`` ``#> tau_sqr[6,5] 7.550537e-04 0.003786986`` ``#> tau_sqr[6,6] 1.102459e-02 0.031718209`` ``#> i_sqr[1,1] 9.975182e-01 0.998826356`` ``#> i_sqr[2,1] 9.970524e-01 0.998508430`` ``#> i_sqr[3,1] 8.999442e-01 0.958593790`` ``#> i_sqr[4,1] 8.531073e-01 0.922697432`` ``#> i_sqr[5,1] 4.580351e-01 0.736897016`` ``#> i_sqr[6,1] 9.323350e-01 0.975081138`

- The fixed part of the random-effects model gives pooled means
  $`\boldsymbol{\alpha} = \mathbb{E} \left[ \mathrm{Vec} \left( \boldsymbol{\mu}, \boldsymbol{\beta} \right)  \right]`$.
- The random part yields between-person covariances
  ($`\boldsymbol{\tau}^{2}`$) quantifying heterogeneity in set point
  ($`\boldsymbol{\mu}`$) and dynamics ($`\boldsymbol{\beta}`$) across
  individuals.

`means`` ``<-`` `[`extract`](https://github.com/jeksterslab/metaDyn/reference/extract.md)`(``metavar``, what ``=`` ``"alpha"``)`` ``means`` ``#> alpha`` ``#> y1 2.97332579`` ``#> y2 2.18436040`` ``#> y3 0.27203469`` ``#> y4 -0.06592972`` ``#> y5 -0.05446716`` ``#> y6 0.24005611`` ``covariances`` ``<-`` `[`extract`](https://github.com/jeksterslab/metaDyn/reference/extract.md)`(``metavar``, what ``=`` ``"tau_sqr"``)`` ``covariances`` ``#> y1 y2 y3 y4 y5 y6`` ``#> y1 1.2395561 0.5388941 0.000000000 0.0000000000 0.0000000000 0.000000000`` ``#> y2 0.5388941 1.1039513 0.000000000 0.0000000000 0.0000000000 0.000000000`` ``#> y3 0.0000000 0.0000000 0.013800790 0.0020642212 0.0013428925 0.003728829`` ``#> y4 0.0000000 0.0000000 0.002064221 0.0098902717 0.0007375405 0.005328959`` ``#> y5 0.0000000 0.0000000 0.001342893 0.0007375405 0.0012891901 0.002271020`` ``#> y6 0.0000000 0.0000000 0.003728829 0.0053289595 0.0022710199 0.021371401`

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
