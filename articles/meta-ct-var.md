# Meta-Analysis of Continuous-Time Vector Autoregressive Model Estimates

## Dynamics Description

The *Stable Reciprocal Regulation* process is modeled as a
**continuous-time bivariate dynamic system** in which two latent
psychological constructs (e.g., negative and positive affect) mutually
influence each other through their instantaneous rates of change. The
dynamics are governed by the continuous-time **drift matrix** whose
**negative diagonal elements** indicate substantial self-regulation:
deviations from a person’s equilibrium are pulled back toward baseline
at a moderate-to-strong rate. The **negative off-diagonal elements**
reflect *reciprocal inhibition*: higher levels of one construct are
associated with an increased instantaneous tendency for the other
construct to decline. Because the cross-effects are modest relative to
the self-regulatory terms, the system exhibits **stable, non-oscillatory
relaxation** back toward person-specific equilibrium points.

Individuals are allowed to differ in both their self-regulatory rates
and the strength of these antagonistic couplings. Stochastic
disturbances enter through the **process-noise covariance**, permitting
small (potentially correlated) departures from the deterministic drift,
while measurement errors are assumed to be minimal and comparable across
indicators. Overall, this pattern captures a psychologically plausible
mechanism of *reciprocal inhibition* in which short-term perturbations
in one component (e.g., negative affect) are naturally counteracted by
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
  \mathrm{d} \boldsymbol{\eta}_{i, t} = \left( \boldsymbol{\alpha}_{i} + \boldsymbol{\beta}_{i} \boldsymbol{\eta}_{i, t} \right) \mathrm{d} t + \boldsymbol{\Psi}_{i}^{\frac{1}{2}} \mathrm{d} \mathbf{W}_{i, t}
\end{equation}
```
where $`\mathrm{d}\boldsymbol{W}`$ is a Wiener process or Brownian
motion, which represents random fluctuations,
$`\boldsymbol{\eta}_{i, t}`$ is a random variable, and
$`\boldsymbol{\beta}_{i}`$, and $`\boldsymbol{\Psi}`$ are model
parameters. Here, $`\boldsymbol{\eta}_{i, t}`$ is a vector of latent
variables at time $`t`$ and individual $`i`$. $`\boldsymbol{\beta}_{i}`$
is a matrix of auto and cross effect coefficients for individual $`i`$,
and $`\boldsymbol{\Psi}`$ the covariance matrix of the process noise.

### Alternative Parameterization

An alternative parameterization of the dynamic structure that directly
estimates the set-point vector $`\boldsymbol{\mu}_{i}`$ is given by
``` math
\begin{equation}
  \mathrm{d} \boldsymbol{\eta}_{i, t} = \boldsymbol{\beta}_{i} \left( \boldsymbol{\eta}_{i, t} - \boldsymbol{\mu}_{i} \right) \mathrm{d} t + \boldsymbol{\Psi}_{i}^{\frac{1}{2}} \mathrm{d} \mathbf{W}_{i, t} .
\end{equation}
```

Algebraic manipulation of the equation results in the following
``` math
\begin{equation}
  \mathrm{d} \boldsymbol{\eta}_{i, t} = \boldsymbol{\beta}_{i} \left( \boldsymbol{\eta}_{i, t} - \boldsymbol{\mu}_{i} \right) \mathrm{d} t + \boldsymbol{\Psi}_{i}^{\frac{1}{2}} \mathrm{d} \mathbf{W}_{i, t} ,
\end{equation}
```
where we can see that the intercept vector $`\boldsymbol{\alpha}_{i}`$
is implied by $`- \boldsymbol{\beta}_{i} \boldsymbol{\mu}_{i}`$.

## Data Generation

### Notation

Let $`t = 1000`$ be the number of time points and $`n = 100`$ be the
number of individuals. We simulate a total of time $`= 1000`$ points per
individual.

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
generate random set point vectors and transition matrices from the
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
#> [1,] 0.5037766 0.1478276
#> [2,] 0.1478276 0.5694164
# first sigma0_l in the list of length n
sigma0_l[[1]] # sigma0_l <- t(chol(sigma0))
#>           [,1]      [,2]
#> [1,] 0.7097722 0.0000000
#> [2,] 0.2082747 0.7252848
# first alpha in the list of length n
alpha[[1]]
#> [1] 3.312112 3.842713
# first beta in the list of length n
beta[[1]]
#>           [,1]      [,2]
#> [1,] -1.235608 -0.186227
#> [2,] -0.169327 -1.325864
# first psi in the list of length n
psi[[1]]
#>      [,1] [,2]
#> [1,] 1.30 0.57
#> [2,] 0.57 1.56
psi_l[[1]] # psi_l <- t(chol(psi))
#>           [,1]     [,2]
#> [1,] 1.1401754 0.000000
#> [2,] 0.4999231 1.144586
# first mu (set-point) in the list of length n
mu[[1]]
#> [1] 2.287769 2.606098
```

### Visualizing the Dynamics Without Process Noise and Measurement Error (n = 5 with Different Initial Condition)

![](fig-vignettes-meta-ct-var-no-error-1.png)![](fig-vignettes-meta-ct-var-no-error-2.png)

### Using the `SimSSMLinSDEFixed` Function from the `simStateSpace` Package to Simulate Data

> **Note:** The `SimSSMLinSDEFixed` function uses a different set of
> parameter names. See
> [`help(SimSSMLinSDEFixed)`](https://github.com/jeksterslab/simStateSpace/reference/SimSSMLinSDEFixed.html)
> for more details.

``` r

library(simStateSpace)
sim <- SimSSMLinSDEIVary(
  n = n,
  time = time,
  delta_t = 0.1,
  mu0 = mu0,
  sigma0_l = sigma0_l,
  iota = alpha,
  phi = beta,
  sigma_l = psi_l,
  nu = nu,
  lambda = lambda,
  theta_l = theta_l
)
data <- as.data.frame(sim)
head(data)
#>   id time       y1       y2
#> 1  1  0.0 2.523315 3.353060
#> 2  1  0.1 2.309397 3.020899
#> 3  1  0.2 2.600318 3.078564
#> 4  1  0.3 2.698363 2.929001
#> 5  1  0.4 2.444623 2.861375
#> 6  1  0.5 3.190602 2.891174
summary(data)
#>        id              time             y1               y2        
#>  Min.   :  1.00   Min.   : 0.00   Min.   :-2.335   Min.   :-1.951  
#>  1st Qu.: 25.75   1st Qu.:24.98   1st Qu.: 2.093   1st Qu.: 1.261  
#>  Median : 50.50   Median :49.95   Median : 2.990   Median : 2.187  
#>  Mean   : 50.50   Mean   :49.95   Mean   : 2.975   Mean   : 2.185  
#>  3rd Qu.: 75.25   3rd Qu.:74.92   3rd Qu.: 3.873   3rd Qu.: 3.081  
#>  Max.   :100.00   Max.   :99.90   Max.   : 8.893   Max.   : 6.970
plot(sim)
```

![](fig-vignettes-meta-ct-var-error-1.png)![](fig-vignettes-meta-ct-var-error-2.png)

## Stage 1: Person-Specific VAR Model

``` r

library(OpenMx)
library(fitVARMxID)
```

The `FitVARMxID` function fits a VAR model on each individual $`i`$.

``` r

fit <- FitVARMxID(
  data = data,
  observed = c("y1", "y2"),
  id = "id",
  ct = TRUE,
  time = "time",
  center = TRUE
)
```

``` r

summary(
  fit,
  means = TRUE
)
#> Call:
#> FitVARMxID(data = data, observed = c("y1", "y2"), id = "id", 
#>     time = "time", ct = TRUE, center = TRUE)
#> 
#> Convergence:
#> 100.0%
#> 
#> Means of the estimated paramaters per individual.
#>   mu[1,1]   mu[2,1] beta[1,1] beta[2,1] beta[1,2] beta[2,2]  psi[1,1]  psi[2,1] 
#>    2.9743    2.1847   -1.3528   -0.1673   -0.1959   -1.4373    1.3040    0.5685 
#>  psi[2,2] 
#>    1.5566
```

## Stage 2: Random-Effects Meta-Analysis of Person-Specific Dynamics and Means

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
library(metaDyn)
metavar <- MetaVARMx(
  object = fit,
  cov_dyn = FALSE,
  tau_sqr_l_free = tau_sqr_l_free,
  tau_sqr_l_values = tau_sqr_l_values
)
```

``` r

summary(metavar)
#> Call:
#> MetaVARMx(object = fit, tau_sqr_l_free = tau_sqr_l_free, tau_sqr_l_values = tau_sqr_l_values, 
#>     cov_dyn = FALSE)
#> 
#> Status code:
#> 0
#> 
#> CI type:
#> "normal"
#> 
#>                  est     se         z      p    2.5%   97.5%
#> alpha[1,1]    2.9744 0.1120   26.5604 0.0000  2.7549  3.1939
#> alpha[2,1]    2.1851 0.1061   20.5982 0.0000  1.9772  2.3930
#> alpha[3,1]   -1.3135 0.0239  -54.9339 0.0000 -1.3604 -1.2666
#> alpha[4,1]   -0.1676 0.0227   -7.3936 0.0000 -0.2120 -0.1231
#> alpha[5,1]   -0.1945 0.0187  -10.3817 0.0000 -0.2312 -0.1578
#> alpha[6,1]   -1.3984 0.0252  -55.4547 0.0000 -1.4478 -1.3489
#> tau_sqr[1,1]  1.2467 0.1775    7.0247 0.0000  0.8989  1.5946
#> tau_sqr[2,1]  0.5415 0.1285    4.2133 0.0000  0.2896  0.7934
#> tau_sqr[2,2]  1.1174 0.1577    7.0845 0.0000  0.8083  1.4265
#> tau_sqr[3,3]  0.0216 0.0071    3.0299 0.0024  0.0076  0.0356
#> tau_sqr[4,3]  0.0017 0.0047    0.3587 0.7198 -0.0075  0.0108
#> tau_sqr[5,3]  0.0001 0.0040    0.0316 0.9748 -0.0078  0.0080
#> tau_sqr[6,3] -0.0050 0.0058   -0.8681 0.3853 -0.0163  0.0063
#> tau_sqr[4,4]  0.0096 0.0059    1.6209 0.1050 -0.0020  0.0212
#> tau_sqr[5,4]  0.0039 0.0030    1.3023 0.1928 -0.0020  0.0098
#> tau_sqr[6,4] -0.0046 0.0054   -0.8424 0.3995 -0.0152  0.0060
#> tau_sqr[5,5]  0.0041 0.0035    1.1897 0.2342 -0.0027  0.0109
#> tau_sqr[6,5]  0.0039 0.0041    0.9542 0.3400 -0.0041  0.0118
#> tau_sqr[6,6]  0.0259 0.0078    3.3250 0.0009  0.0106  0.0412
#> i_sqr[1,1]    0.9950 0.0007 1397.1450 0.0000  0.9936  0.9964
#> i_sqr[2,1]    0.9932 0.0010 1038.1662 0.0000  0.9913  0.9951
#> i_sqr[3,1]    0.3951 0.0787    5.0232 0.0000  0.2409  0.5493
#> i_sqr[4,1]    0.2663 0.0999    2.6661 0.0077  0.0705  0.4621
#> i_sqr[5,1]    0.1684 0.0960    1.7552 0.0792 -0.0197  0.3565
#> i_sqr[6,1]    0.4450 0.0723    6.1563 0.0000  0.3034  0.5867
```

### Normal Theory Confidence Intervals

``` r

confint(metavar, level = 0.95)
#>                     2.5 %       97.5 %
#> alpha[1,1]    2.754884057  3.193857583
#> alpha[2,1]    1.977166171  2.392997102
#> alpha[3,1]   -1.360361699 -1.266634180
#> alpha[4,1]   -0.211971209 -0.123137203
#> alpha[5,1]   -0.231196553 -0.157764459
#> alpha[6,1]   -1.447791742 -1.348945259
#> tau_sqr[1,1]  0.898881136  1.594579410
#> tau_sqr[2,1]  0.289612006  0.793425685
#> tau_sqr[2,2]  0.808257178  1.426525394
#> tau_sqr[3,3]  0.007629055  0.035580870
#> tau_sqr[4,3] -0.007465183  0.010809554
#> tau_sqr[5,3] -0.007793433  0.008049002
#> tau_sqr[6,3] -0.016283531  0.006286548
#> tau_sqr[4,4] -0.002009185  0.021220691
#> tau_sqr[5,4] -0.001983889  0.009841723
#> tau_sqr[6,4] -0.015168121  0.006048596
#> tau_sqr[5,5] -0.002677109  0.010946654
#> tau_sqr[6,5] -0.004084324  0.011833761
#> tau_sqr[6,6]  0.010649363  0.041230450
#> i_sqr[1,1]    0.993576205  0.996367768
#> i_sqr[2,1]    0.991335776  0.995085961
#> i_sqr[3,1]    0.240935212  0.549253310
#> i_sqr[4,1]    0.070527851  0.462061436
#> i_sqr[5,1]   -0.019652611  0.356545330
#> i_sqr[6,1]    0.303358017  0.586733849
```

``` r

confint(metavar, level = 0.99)
#>                     0.5 %       99.5 %
#> alpha[1,1]    2.685916317  3.262825322
#> alpha[2,1]    1.911834397  2.458328876
#> alpha[3,1]   -1.375087360 -1.251908519
#> alpha[4,1]   -0.225928042 -0.109180369
#> alpha[5,1]   -0.242733571 -0.146227441
#> alpha[6,1]   -1.463321649 -1.333415351
#> tau_sqr[1,1]  0.789579016  1.703881530
#> tau_sqr[2,1]  0.210457140  0.872580551
#> tau_sqr[2,2]  0.711120201  1.523662370
#> tau_sqr[3,3]  0.003237507  0.039972418
#> tau_sqr[4,3] -0.010336353  0.013680724
#> tau_sqr[5,3] -0.010282460  0.010538029
#> tau_sqr[6,3] -0.019829548  0.009832564
#> tau_sqr[4,4] -0.005658863  0.024870369
#> tau_sqr[5,4] -0.003841828  0.011699661
#> tau_sqr[6,4] -0.018501509  0.009381984
#> tau_sqr[5,5] -0.004817557  0.013087102
#> tau_sqr[6,5] -0.006585237  0.014334674
#> tau_sqr[6,6]  0.005844727  0.046035087
#> i_sqr[1,1]    0.993137619  0.996806355
#> i_sqr[2,1]    0.990746579  0.995675157
#> i_sqr[3,1]    0.192494928  0.597693594
#> i_sqr[4,1]    0.009013466  0.523575821
#> i_sqr[5,1]   -0.078757592  0.415650311
#> i_sqr[6,1]    0.258836447  0.631255419
```

### Robust Confidence Intervals

``` r

confint(metavar, level = 0.95, robust = TRUE)
#>                      2.5 %       97.5 %
#> alpha[1,1]    2.7537489604  3.194992679
#> alpha[2,1]    1.9759518792  2.394211394
#> alpha[3,1]   -1.3588674820 -1.268128397
#> alpha[4,1]   -0.2123668782 -0.122741533
#> alpha[5,1]   -0.2310729994 -0.157888013
#> alpha[6,1]   -1.4453893048 -1.351347696
#> tau_sqr[1,1]  0.9084167397  1.585043807
#> tau_sqr[2,1]  0.2934569960  0.789580695
#> tau_sqr[2,2]  0.8505421607  1.384240411
#> tau_sqr[3,3]  0.0060477632  0.037162162
#> tau_sqr[4,3] -0.0069557177  0.010300088
#> tau_sqr[5,3] -0.0076598547  0.007915424
#> tau_sqr[6,3] -0.0153198902  0.005322906
#> tau_sqr[4,4] -0.0003311056  0.019542611
#> tau_sqr[5,4] -0.0010859415  0.008943775
#> tau_sqr[6,4] -0.0132340711  0.004114546
#> tau_sqr[5,5] -0.0017118526  0.009981397
#> tau_sqr[6,5] -0.0031210890  0.010870526
#> tau_sqr[6,6]  0.0050732355  0.046806578
#> i_sqr[1,1]    0.9936144676  0.996329506
#> i_sqr[2,1]    0.9915998882  0.994821848
#> i_sqr[3,1]    0.2234487839  0.566739739
#> i_sqr[4,1]    0.0900913529  0.442497934
#> i_sqr[5,1]    0.0122705694  0.324622150
#> i_sqr[6,1]    0.2267817671  0.663310098
```

``` r

confint(metavar, level = 0.99, robust = TRUE)
#>                     0.5 %       99.5 %
#> alpha[1,1]    2.684424548  3.264317091
#> alpha[2,1]    1.910238547  2.459924726
#> alpha[3,1]   -1.373123626 -1.253872254
#> alpha[4,1]   -0.226448040 -0.108660371
#> alpha[5,1]   -0.242571194 -0.146389818
#> alpha[6,1]   -1.460164312 -1.336572688
#> tau_sqr[1,1]  0.802110923  1.691349624
#> tau_sqr[2,1]  0.215510314  0.867527377
#> tau_sqr[2,2]  0.766692089  1.468090483
#> tau_sqr[3,3]  0.001159337  0.042050588
#> tau_sqr[4,3] -0.009666801  0.013011172
#> tau_sqr[5,3] -0.010106908  0.010362477
#> tau_sqr[6,3] -0.018563109  0.008566125
#> tau_sqr[4,4] -0.003453493  0.022664999
#> tau_sqr[5,4] -0.002661724  0.010519558
#> tau_sqr[6,4] -0.015959736  0.006840212
#> tau_sqr[5,5] -0.003548995  0.011818540
#> tau_sqr[6,5] -0.005319331  0.013068768
#> tau_sqr[6,6] -0.001483548  0.053363361
#> i_sqr[1,1]    0.993187904  0.996756069
#> i_sqr[2,1]    0.991093682  0.995328055
#> i_sqr[3,1]    0.169513866  0.620674657
#> i_sqr[4,1]    0.034724266  0.497865021
#> i_sqr[5,1]   -0.036803421  0.373696139
#> i_sqr[6,1]    0.158198196  0.731893670
```

- The fixed part of the random-effects model gives pooled means
  $`\boldsymbol{\alpha} = \mathbb{E} \left[ \mathrm{Vec} \left( \boldsymbol{\mu}, \boldsymbol{\beta} \right)  \right]`$.
- The random part yields between-person covariances
  ($`\boldsymbol{\tau}^{2}`$) quantifying heterogeneity in set-point
  ($`\boldsymbol{\mu}`$) and dynamics ($`\boldsymbol{\beta}`$) across
  individuals.

``` r

means <- extract(metavar, what = "alpha")
means
#>         alpha
#> y1  2.9743708
#> y2  2.1850816
#> y3 -1.3134979
#> y4 -0.1675542
#> y5 -0.1944805
#> y6 -1.3983685
covariances <- extract(metavar, what = "tau_sqr")
covariances
#>            y1           y2            y3           y4           y5           y6
#> y1 1.24673027  0.541518845  0.0255140119  0.018167265 0.0347879308  0.029621672
#> y2 0.54151885  1.117391286 -0.0082835061 -0.033976084 0.0242975818  0.036579585
#> y3 0.02551401 -0.008283506  0.0216049626  0.001672185 0.0001277845 -0.004998492
#> y4 0.01816726 -0.033976084  0.0016721854  0.009605753 0.0039289167 -0.004559762
#> y5 0.03478793  0.024297582  0.0001277845  0.003928917 0.0041347723  0.003874719
#> y6 0.02962167  0.036579585 -0.0049984919 -0.004559762 0.0038747186  0.025939907
```

## Summary

This vignette demonstrates a two-stage hierarchical estimation approach
for dynamic systems: 1. individual-level CT-VAR estimation, and\
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
