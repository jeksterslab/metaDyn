.CIMCMeta <- function(object,
                      alpha,
                      nrep = 20000L,
                      seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  model <- object$output

  mat_names <- c(
    names(model$algebras),
    names(model$matrices)
  )

  .EvalByNameVec <- function(name,
                             model) {
    out <- OpenMx::mxEvalByName(
      name = name,
      model = model,
      compute = TRUE
    )

    out_vec <- c(out)

    if (grepl(pattern = "^ie_", x = name) && length(out_vec) == 1L) {
      out_names <- name
    } else {
      out_names <- rownames(out)

      if (is.null(out_names) || length(out_names) != length(out_vec)) {
        out_names <- names(out_vec)
      }

      if (is.null(out_names) || anyNA(out_names) || any(!nzchar(out_names))) {
        out_names <- paste0(
          name,
          "[",
          seq_along(out_vec),
          "]"
        )
      }
    }

    names(out_vec) <- out_names
    out_vec
  }

  .Rmvn <- function(n,
                    mu,
                    sigma) {
    sigma <- 0.5 * (
      sigma + t(sigma)
    )

    eig <- eigen(
      x = sigma,
      symmetric = TRUE
    )

    eig$values <- ifelse(
      test = eig$values < 0,
      yes = 0,
      no = eig$values
    )

    root <- eig$vectors %*% diag(
      x = sqrt(eig$values),
      nrow = length(eig$values),
      ncol = length(eig$values)
    )

    z <- matrix(
      data = stats::rnorm(
        n = n * length(mu)
      ),
      nrow = n,
      ncol = length(mu)
    )

    out <- sweep(
      x = z %*% t(root),
      MARGIN = 2,
      STATS = mu,
      FUN = "+"
    )

    colnames(out) <- names(mu)
    out
  }

  .CIMC <- function(est,
                    draws,
                    alpha,
                    theta = 0) {
    probs <- .ProbsofAlpha(
      alpha = alpha
    )

    se <- apply(
      X = draws,
      MARGIN = 2,
      FUN = stats::sd,
      na.rm = TRUE
    )

    stat <- (
      est - theta
    ) / se

    p_lower <- colMeans(
      x = sweep(
        x = draws,
        MARGIN = 2,
        STATS = theta,
        FUN = "<="
      ),
      na.rm = TRUE
    )

    p_upper <- colMeans(
      x = sweep(
        x = draws,
        MARGIN = 2,
        STATS = theta,
        FUN = ">="
      ),
      na.rm = TRUE
    )

    p <- 2 * pmin(
      p_lower,
      p_upper
    )

    p <- pmin(
      p,
      1
    )

    ci <- apply(
      X = draws,
      MARGIN = 2,
      FUN = stats::quantile,
      probs = probs,
      na.rm = TRUE,
      names = FALSE
    )

    ci <- t(ci)

    out <- cbind(
      est = est,
      se = se,
      z = stat,
      p = p,
      ci
    )

    colnames(out) <- c(
      "est",
      "se",
      "z",
      "p",
      paste0(
        probs * 100,
        "%"
      )
    )

    rownames(out) <- names(est)

    out
  }

  target_names <- c(
    "alpha_vec",
    "gamma_vec",
    "mu_x_vec",
    "beta_vec",
    "kappa_vec",
    "phi_vec",
    "omega_vec",
    "psi_vec",
    "sigma_x_vec",
    "tau_sqr_vec",
    "i_sqr_vec",
    "direct_vec",
    "indirect_vec",
    "total_vec"
  )

  ie_xyz <- grep(
    pattern = "^ie_x\\d+_y\\d+_z\\d+$",
    x = mat_names,
    value = TRUE
  )

  target_names <- c(
    target_names,
    ie_xyz
  )

  target_names <- target_names[
    target_names %in% mat_names
  ]

  est_list <- lapply(
    X = target_names,
    FUN = .EvalByNameVec,
    model = model
  )

  names(est_list) <- target_names

  pars <- stats::coef(
    object = model
  )

  vcov_pars <- stats::vcov(
    object = model
  )

  if (is.null(names(pars)) || any(!nzchar(names(pars)))) {
    stop(
      "Free parameters must have names for Monte Carlo confidence intervals.",
      call. = FALSE
    )
  }

  if (is.null(rownames(vcov_pars)) || is.null(colnames(vcov_pars))) {
    rownames(vcov_pars) <- colnames(vcov_pars) <- names(pars)
  }

  missing_pars <- setdiff(
    x = names(pars),
    y = rownames(vcov_pars)
  )

  if (length(missing_pars) > 0) {
    stop(
      paste(
        "The parameter covariance matrix is missing the following parameters:",
        paste(
          missing_pars,
          collapse = ", "
        )
      ),
      call. = FALSE
    )
  }

  vcov_pars <- vcov_pars[
    names(pars),
    names(pars),
    drop = FALSE
  ]

  if (any(!is.finite(vcov_pars))) {
    stop(
      "The parameter covariance matrix contains non-finite values.",
      call. = FALSE
    )
  }

  vcov_pars <- 0.5 * (
    vcov_pars + t(vcov_pars)
  )

  par_draws <- .Rmvn(
    n = nrep,
    mu = pars,
    sigma = vcov_pars
  )

  draws_list <- lapply(
    X = est_list,
    FUN = function(x) {
      out <- matrix(
        data = NA_real_,
        nrow = nrep,
        ncol = length(x)
      )
      colnames(out) <- names(x)
      out
    }
  )

  ok <- rep(
    x = TRUE,
    times = nrep
  )

  for (i in seq_len(nrep)) {
    values_i <- as.numeric(
      par_draws[i, ]
    )
    names(values_i) <- names(pars)

    model_i <- tryCatch(
      expr = OpenMx::omxSetParameters(
        model = model,
        labels = names(values_i),
        values = values_i
      ),
      error = function(e) {
        NULL
      }
    )

    if (is.null(model_i)) {
      ok[i] <- FALSE
    } else {
      for (j in seq_along(target_names)) {
        val <- tryCatch(
          expr = .EvalByNameVec(
            name = target_names[j],
            model = model_i
          ),
          error = function(e) {
            rep(
              x = NA_real_,
              times = length(est_list[[j]])
            )
          }
        )

        draws_list[[j]][i, ] <- val
      }
    }
  }
  # -------------------------------------------------------------

  draws_list <- lapply(
    X = draws_list,
    FUN = function(x) {
      x[
        ok, ,
        drop = FALSE
      ]
    }
  )

  if (!any(ok)) {
    stop(
      "No valid Monte Carlo draws were generated.",
      call. = FALSE
    )
  }

  out <- list()

  if ("alpha_vec" %in% target_names) {
    out$y0 <- .CIMC(
      est = est_list$alpha_vec,
      draws = draws_list$alpha_vec,
      alpha = alpha
    )
  }

  if ("gamma_vec" %in% target_names) {
    out$y1 <- .CIMC(
      est = est_list$gamma_vec,
      draws = draws_list$gamma_vec,
      alpha = alpha
    )
  }

  if ("mu_x_vec" %in% target_names) {
    out$x0 <- .CIMC(
      est = est_list$mu_x_vec,
      draws = draws_list$mu_x_vec,
      alpha = alpha
    )
  }

  if ("beta_vec" %in% target_names) {
    out$yy <- .CIMC(
      est = est_list$beta_vec,
      draws = draws_list$beta_vec,
      alpha = alpha
    )
  }

  if ("kappa_vec" %in% target_names) {
    out$z0 <- .CIMC(
      est = est_list$kappa_vec,
      draws = draws_list$kappa_vec,
      alpha = alpha
    )
  }

  if ("phi_vec" %in% target_names) {
    out$z1 <- .CIMC(
      est = est_list$phi_vec,
      draws = draws_list$phi_vec,
      alpha = alpha
    )
  }

  if ("omega_vec" %in% target_names) {
    out$zx <- .CIMC(
      est = est_list$omega_vec,
      draws = draws_list$omega_vec,
      alpha = alpha
    )
  }

  if ("psi_vec" %in% target_names) {
    out$psi <- .CIMC(
      est = est_list$psi_vec,
      draws = draws_list$psi_vec,
      alpha = alpha
    )
  }

  if ("sigma_x_vec" %in% target_names) {
    out$sx <- .CIMC(
      est = est_list$sigma_x_vec,
      draws = draws_list$sigma_x_vec,
      alpha = alpha
    )
  }

  if ("tau_sqr_vec" %in% target_names) {
    out$t2 <- .CIMC(
      est = est_list$tau_sqr_vec,
      draws = draws_list$tau_sqr_vec,
      alpha = alpha
    )
  }

  if ("i_sqr_vec" %in% target_names) {
    out$i2 <- .CIMC(
      est = est_list$i_sqr_vec,
      draws = draws_list$i_sqr_vec,
      alpha = alpha
    )
  }

  if ("direct_vec" %in% target_names) {
    out$direct_vec <- .CIMC(
      est = est_list$direct_vec,
      draws = draws_list$direct_vec,
      alpha = alpha
    )
  }

  if ("indirect_vec" %in% target_names) {
    out$indirect_vec <- .CIMC(
      est = est_list$indirect_vec,
      draws = draws_list$indirect_vec,
      alpha = alpha
    )
  }

  if ("total_vec" %in% target_names) {
    out$total_vec <- .CIMC(
      est = est_list$total_vec,
      draws = draws_list$total_vec,
      alpha = alpha
    )
  }

  if (length(ie_xyz) > 0) {
    out$ie_xyz <- do.call(
      what = "rbind",
      args = lapply(
        X = ie_xyz,
        FUN = function(name) {
          .CIMC(
            est = est_list[[name]],
            draws = draws_list[[name]],
            alpha = alpha
          )
        }
      )
    )
  }

  out <- out[
    !sapply(
      X = out,
      FUN = is.null
    )
  ]

  do.call(
    what = "rbind",
    args = out
  )
}
