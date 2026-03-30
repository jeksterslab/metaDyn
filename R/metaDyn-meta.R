#' Fit Multivariate Meta-Analysis
#'
#' This function estimates fixed-, random-, or mixed-effects
#' meta-analytic parameters using per-individual coefficient estimates
#' and their sampling variance-covariance matrices.
#' Optionally, it fits distal-outcome models
#' in which between-person outcomes are regressed on
#' between-person covariates and the meta-analyzed parameters/effect sizes.
#'
#' @author Ivan Jacob Agaloos Pesigan
#'
#' @param y A list.
#'   Each element of the list is a numeric vector
#'   of estimated coefficients.
#' @param v A list.
#'   Each element of the list
#'   is a sampling variance-covariance matrix of `y`.
#' @param x An optional list.
#'   Each element of the list is a numeric vector of covariates.
#' @param z An optional list.
#'   Each element of the list is a numeric vector of distal outcomes.
#' @param random Logical.
#'   If `random = TRUE`,
#'   estimates random effects.
#'   If `random = FALSE`,
#'   `tau_sqr` is a null matrix.
#' @param alpha_free Logical vector.
#'   Optional vector of free (`TRUE`) parameters for `alpha`.
#' @param alpha_values Numeric vector.
#'   Optional vector of starting values for `alpha`.
#' @param alpha_lbound Numeric vector.
#'   Optional vector of lower bound values for `alpha`.
#' @param alpha_ubound Numeric vector.
#'   Optional vector of upper bound values for `alpha`.
#' @param tau_sqr_diag Logical.
#'   If `tau_sqr_diag = TRUE`,
#'   `tau_sqr` is a diagonal matrix.
#'   If `tau_sqr_diag = FALSE`,
#'   `tau_sqr` is a symmetric matrix.
#' @param tau_sqr_d_free Logical vector
#'   indicating free/fixed status of the elements of `tau_sqr_d`.
#'   If `NULL`, all element of `tau_sqr_d` are free.
#' @param tau_sqr_d_values Numeric vector
#'   with starting values for `tau_sqr_d`.
#'   If `NULL`, defaults to a vector of ones.
#' @param tau_sqr_d_lbound Numeric vector
#'   with lower bounds for `tau_sqr_d`.
#'   If `NULL`, no lower bounds are set.
#' @param tau_sqr_d_ubound Numeric vector
#'   with upper bounds for `tau_sqr_d`.
#'   If `NULL`, no upper bounds are set.
#' @param tau_sqr_l_free Logical matrix
#'   indicating which strictly-lower-triangular elements
#'   of `tau_sqr_l` are free.
#'   Ignored if `tau_sqr_diag = TRUE`.
#' @param tau_sqr_l_values Numeric matrix
#'   of starting values for the strictly-lower-triangular elements
#'   of `tau_sqr_l`.
#'   If `NULL`, defaults to a null matrix.
#' @param tau_sqr_l_lbound Numeric matrix
#'   with lower bounds for `tau_sqr_l`.
#'   If `NULL`, no lower bounds are set.
#' @param tau_sqr_l_ubound Numeric matrix
#'   with upper bounds for `tau_sqr_l`.
#'   If `NULL`, no upper bounds are set.
#' @param gamma_free Logical matrix.
#'   Optional matrix of free (`TRUE`) parameters for `gamma`.
#' @param gamma_values Numeric matrix.
#'   Optional matrix of starting values for `gamma`.
#' @param gamma_lbound Numeric matrix.
#'   Optional matrix of lower bound values for `gamma`.
#' @param gamma_ubound Numeric matrix.
#'   Optional matrix of upper bound values for `gamma`.
#' @param kappa_free Logical vector.
#'   Optional vector of free (`TRUE`) parameters for `kappa`.
#' @param kappa_values Numeric vector.
#'   Optional vector of starting values for `kappa`.
#' @param kappa_lbound Numeric vector.
#'   Optional vector of lower bound values for `kappa`.
#' @param kappa_ubound Numeric vector.
#'   Optional vector of upper bound values for `kappa`.
#' @param phi_free Logical matrix.
#'   Optional matrix of free (`TRUE`) parameters for `phi`.
#' @param phi_values Numeric matrix.
#'   Optional matrix of starting values for `phi`.
#' @param phi_lbound Numeric matrix.
#'   Optional matrix of lower bound values for `phi`.
#' @param phi_ubound Numeric matrix.
#'   Optional matrix of upper bound values for `phi`.
#' @param omega_free Logical matrix.
#'   Optional matrix of free (`TRUE`) parameters for `omega`.
#' @param omega_values Numeric matrix.
#'   Optional matrix of starting values for `omega`.
#' @param omega_lbound Numeric matrix.
#'   Optional matrix of lower bound values for `omega`.
#' @param omega_ubound Numeric matrix.
#'   Optional matrix of upper bound values for `omega`.
#' @param psi_diag Logical.
#'   If `psi_diag = TRUE`,
#'   `psi` is a diagonal matrix.
#'   If `psi_diag = FALSE`,
#'   `psi` is a symmetric matrix.
#' @param psi_d_free Logical vector
#'   indicating free/fixed status of the elements of `psi_d`.
#'   If `NULL`, all element of `psi_d` are free.
#' @param psi_d_values Numeric vector
#'   with starting values for `psi_d`.
#'   If `NULL`, defaults to a vector of ones.
#' @param psi_d_lbound Numeric vector
#'   with lower bounds for `psi_d`.
#'   If `NULL`, no lower bounds are set.
#' @param psi_d_ubound Numeric vector
#'   with upper bounds for `psi_d`.
#'   If `NULL`, no upper bounds are set.
#' @param psi_l_free Logical matrix
#'   indicating which strictly-lower-triangular elements
#'   of `psi_l` are free.
#'   Ignored if `psi_diag = TRUE`.
#' @param psi_l_values Numeric matrix
#'   of starting values for the strictly-lower-triangular elements
#'   of `psi_l`.
#'   If `NULL`, defaults to a null matrix.
#' @param psi_l_lbound Numeric matrix
#'   with lower bounds for `psi_l`.
#'   If `NULL`, no lower bounds are set.
#' @param psi_l_ubound Numeric matrix
#'   with upper bounds for `psi_l`.
#'   If `NULL`, no upper bounds are set.
#' @param check_estimates Logical.
#'   Check elements of `v` for positive definiteness.
#'   If the test fails, the function generates a near positive definite matrix
#'   to replace the original using [Matrix::nearPD()].
#' @param robust Logical.
#'   If `TRUE`,
#'   calculate robust (sandwich) sampling variance-covariance matrix
#'   in stage 2.
#' @param alpha NUmeric.
#'   Alpha for test of significance and confidence intervals.
#' @param tries_explore Integer.
#'   Number of extra tries for the wide exploration phase.
#' @param tries_local Integer.
#'   Number of extra tries for local polishing.
#' @param max_attempts Integer.
#'   Maximum number of remediation attempts
#'   after the first Hessian computation fails the criteria.
#' @param silent Logical.
#'   If `TRUE`, suppresses messages during the model fitting stage.
#' @param seed Random seed for reproducibility.
#' @param ncores Positive integer.
#'   Number of cores to use.
#'
#' @return Returns an object of class `metadynmeta` which is
#'   a list with the following elements:
#'   \describe{
#'     \item{call}{Function call.}
#'     \item{args}{List of function arguments.}
#'     \item{fun}{Function used ("Meta").}
#'     \item{output}{A fitted OpenMx model.}
#'     \item{robust}{Output from [OpenMx::imxRobustSE()]
#'         with argument `details = TRUE` if `robust = TRUE`.}
#'   }
#'
#' @examples
#' \donttest{
#' if (requireNamespace("simStateSpace")) {
#'   # Generate data using the simStateSpace package-------------------------
#'   library(simStateSpace)
#'   set.seed(42)
#'   n <- 5
#'   time <- 100
#'   p <- 2
#'   alpha <- rep(x = 0, times = p)
#'   beta <- 0.50 * diag(p)
#'   psi <- 0.001 * diag(p)
#'   psi_l <- t(chol(psi))
#'   mu0 <- SSMMeanEta(
#'     beta = beta,
#'     alpha = alpha
#'   )
#'   sigma0 <- SSMCovEta(
#'     beta = beta,
#'     psi = psi
#'   )
#'   sigma0_l <- t(chol(sigma0))
#'   sim <- SimSSMVARFixed(
#'     n = n,
#'     time = time,
#'     mu0 = mu0,
#'     sigma0_l = sigma0_l,
#'     alpha = alpha,
#'     beta = beta,
#'     psi_l = psi_l
#'   )
#'   data <- as.data.frame(sim)
#'
#'   # Stage 1---------------------------------------------------------------
#'   library(fitVARMxID)
#'   stage1 <- FitVARMxID(
#'     data = data,
#'     observed = paste0("y", seq_len(p)),
#'     id = "id",
#'     center = TRUE
#'   )
#'   summary(stage1)
#'   # Stage 2---------------------------------------------------------------
#'   # Meta-analyze set point vector and matrix of lagged-effects
#'   y <- coef(
#'     object = stage1,
#'     mu = TRUE,
#'     beta = TRUE,
#'     alpha = FALSE,
#'     nu = FALSE,
#'     psi = FALSE,
#'     theta = FALSE
#'   )
#'   v <- vcov(
#'     object = stage1,
#'     mu = TRUE,
#'     beta = TRUE,
#'     alpha = FALSE,
#'     nu = FALSE,
#'     psi = FALSE,
#'     theta = FALSE
#'   )
#'   library(metaDyn)
#'   stage2 <- Meta(y = y, v = v, random = FALSE)
#'   # Methods for the output of the Meta() function
#'   print(stage2)
#'   summary(stage2)
#'   coef(stage2)
#'   vcov(stage2)
#'   confint(stage2)
#'   extract(stage2, what = "alpha")
#' }
#' }
#'
#' @references
#' Cheung, M. W.-L. (2015).
#' *Meta-analysis: A structural equation modeling approach*.
#' Wiley.
#' \doi{10.1002/9781118957813}
#'
#' Neale, M. C., Hunter, M. D., Pritikin, J. N.,
#' Zahery, M., Brick, T. R., Kirkpatrick, R. M., Estabrook, R.,
#' Bates, T. C., Maes, H. H., & Boker, S. M. (2015).
#' OpenMx 2.0: Extended structural equation and statistical modeling.
#' *Psychometrika*,
#' *81*(2), 535–549.
#' \doi{10.1007/s11336-014-9435-8}
#'
#' @family Meta-Analysis of VAR Functions
#' @keywords metaDyn meta
#' @import OpenMx
#' @importFrom stats coef vcov
#' @export
Meta <- function(y,
                 v,
                 x = NULL,
                 z = NULL,
                 random = TRUE,
                 alpha_free = NULL,
                 alpha_values = NULL,
                 alpha_lbound = NULL,
                 alpha_ubound = NULL,
                 tau_sqr_diag = FALSE,
                 tau_sqr_d_free = NULL,
                 tau_sqr_d_values = NULL,
                 tau_sqr_d_lbound = NULL,
                 tau_sqr_d_ubound = NULL,
                 tau_sqr_l_free = NULL,
                 tau_sqr_l_values = NULL,
                 tau_sqr_l_lbound = NULL,
                 tau_sqr_l_ubound = NULL,
                 gamma_free = NULL,
                 gamma_values = NULL,
                 gamma_lbound = NULL,
                 gamma_ubound = NULL,
                 kappa_free = NULL,
                 kappa_values = NULL,
                 kappa_lbound = NULL,
                 kappa_ubound = NULL,
                 phi_free = NULL,
                 phi_values = NULL,
                 phi_lbound = NULL,
                 phi_ubound = NULL,
                 omega_free = NULL,
                 omega_values = NULL,
                 omega_lbound = NULL,
                 omega_ubound = NULL,
                 psi_diag = FALSE,
                 psi_d_free = NULL,
                 psi_d_values = NULL,
                 psi_d_lbound = NULL,
                 psi_d_ubound = NULL,
                 psi_l_free = NULL,
                 psi_l_values = NULL,
                 psi_l_lbound = NULL,
                 psi_l_ubound = NULL,
                 check_estimates = TRUE,
                 robust = FALSE,
                 alpha = 0.05,
                 seed = NULL,
                 tries_explore = 100,
                 tries_local = 100,
                 max_attempts = 10,
                 silent = FALSE,
                 ncores = NULL) {
  .MetaCheckFree(
    alpha_free = alpha_free,
    alpha_values = alpha_values,
    tau_sqr_diag = tau_sqr_diag,
    tau_sqr_d_free = tau_sqr_d_free,
    tau_sqr_d_values = tau_sqr_d_values,
    tau_sqr_l_free = tau_sqr_l_free,
    tau_sqr_l_values = tau_sqr_l_values,
    gamma_free = gamma_free,
    gamma_values = gamma_values,
    kappa_free = kappa_free,
    kappa_values = kappa_values,
    phi_free = phi_free,
    phi_values = phi_values,
    omega_free = omega_free,
    omega_values = omega_values,
    psi_diag = psi_diag,
    psi_d_free = psi_d_free,
    psi_d_values = psi_d_values,
    psi_l_free = psi_l_free,
    psi_l_values = psi_l_values
  )
  p <- length(y[[1]])
  n <- length(y)
  stopifnot(
    length(v) == length(y)
  )
  if (is.null(x)) {
    m <- NULL
    covariate <- FALSE
  } else {
    m <- length(x[[1]])
    stopifnot(
      length(x) == length(y)
    )
    covariate <- TRUE
  }
  if (is.null(z)) {
    r <- NULL
    distal <- FALSE
  } else {
    r <- length(z[[1]])
    stopifnot(
      length(z) == length(y)
    )
    distal <- TRUE
  }
  if (check_estimates) {
    estimates <- .CheckEstimates(
      y = y,
      v = v,
      p = p,
      ynames = paste0("y", seq_len(p)),
      ncores = ncores
    )
    y <- estimates$y
    v <- estimates$v
  }
  args <- list(
    y = y,
    v = v,
    x = x,
    z = z,
    p = p,
    m = m,
    r = r,
    n = n,
    random = random,
    covariate = covariate,
    distal = distal,
    alpha_free = alpha_free,
    alpha_values = alpha_values,
    alpha_lbound = alpha_lbound,
    alpha_ubound = alpha_ubound,
    tau_sqr_diag = tau_sqr_diag,
    tau_sqr_d_free = tau_sqr_d_free,
    tau_sqr_d_values = tau_sqr_d_values,
    tau_sqr_d_lbound = tau_sqr_d_lbound,
    tau_sqr_d_ubound = tau_sqr_d_ubound,
    tau_sqr_l_free = tau_sqr_l_free,
    tau_sqr_l_values = tau_sqr_l_values,
    tau_sqr_l_lbound = tau_sqr_l_lbound,
    tau_sqr_l_ubound = tau_sqr_l_ubound,
    gamma_free = gamma_free,
    gamma_values = gamma_values,
    gamma_lbound = gamma_lbound,
    gamma_ubound = gamma_ubound,
    kappa_free = kappa_free,
    kappa_values = kappa_values,
    kappa_lbound = kappa_lbound,
    kappa_ubound = kappa_ubound,
    phi_free = phi_free,
    phi_values = phi_values,
    phi_lbound = phi_lbound,
    phi_ubound = phi_ubound,
    omega_free = omega_free,
    omega_values = omega_values,
    omega_lbound = omega_lbound,
    omega_ubound = omega_ubound,
    psi_diag = psi_diag,
    psi_d_free = psi_d_free,
    psi_d_values = psi_d_values,
    psi_d_lbound = psi_d_lbound,
    psi_d_ubound = psi_d_ubound,
    psi_l_free = psi_l_free,
    psi_l_values = psi_l_values,
    psi_l_lbound = psi_l_lbound,
    psi_l_ubound = psi_l_ubound,
    check_estimates = check_estimates,
    robust = robust,
    alpha = alpha,
    seed = seed,
    tries_explore = tries_explore,
    tries_local = tries_local,
    max_attempts = max_attempts,
    silent = silent,
    ncores = ncores
  )
  output <- .MetaFit(
    y = y,
    v = v,
    x = x,
    z = z,
    p = p,
    m = m,
    r = r,
    n = n,
    random = random,
    covariate = covariate,
    distal = distal,
    alpha_free = alpha_free,
    alpha_values = alpha_values,
    alpha_lbound = alpha_lbound,
    alpha_ubound = alpha_ubound,
    tau_sqr_diag = tau_sqr_diag,
    tau_sqr_d_free = tau_sqr_d_free,
    tau_sqr_d_values = tau_sqr_d_values,
    tau_sqr_d_lbound = tau_sqr_d_lbound,
    tau_sqr_d_ubound = tau_sqr_d_ubound,
    tau_sqr_l_free = tau_sqr_l_free,
    tau_sqr_l_values = tau_sqr_l_values,
    tau_sqr_l_lbound = tau_sqr_l_lbound,
    tau_sqr_l_ubound = tau_sqr_l_ubound,
    gamma_free = gamma_free,
    gamma_values = gamma_values,
    gamma_lbound = gamma_lbound,
    gamma_ubound = gamma_ubound,
    kappa_free = kappa_free,
    kappa_values = kappa_values,
    kappa_lbound = kappa_lbound,
    kappa_ubound = kappa_ubound,
    phi_free = phi_free,
    phi_values = phi_values,
    phi_lbound = phi_lbound,
    phi_ubound = phi_ubound,
    omega_free = omega_free,
    omega_values = omega_values,
    omega_lbound = omega_lbound,
    omega_ubound = omega_ubound,
    psi_diag = psi_diag,
    psi_d_free = psi_d_free,
    psi_d_values = psi_d_values,
    psi_d_lbound = psi_d_lbound,
    psi_d_ubound = psi_d_ubound,
    psi_l_free = psi_l_free,
    psi_l_values = psi_l_values,
    psi_l_lbound = psi_l_lbound,
    psi_l_ubound = psi_l_ubound,
    alpha = alpha,
    seed = seed,
    tries_explore = tries_explore,
    tries_local = tries_local,
    max_attempts = max_attempts,
    silent = silent,
    ncores = ncores
  )
  if (robust) {
    utils::capture.output(
      suppressMessages(
        suppressWarnings(
          sandwich <- OpenMx::imxRobustSE(
            model = output,
            details = TRUE
          )
        )
      )
    )
  } else {
    sandwich <- NULL
  }
  out <- list(
    call = match.call(),
    args = args,
    fun = "Meta",
    output = output,
    robust = sandwich
  )
  class(out) <- c(
    "metadynmeta",
    class(out)
  )
  out
}
