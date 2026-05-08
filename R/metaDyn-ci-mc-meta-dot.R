.CIMCMeta <- function(object,
                      alpha,
                      nrep,
                      seed,
                      ncores) {
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
    r <- colSums(
      !is.na(draws)
    )
    se <- apply(
      X = draws,
      MARGIN = 2,
      FUN = function(x) {
        stats::sd(
          x = x,
          na.rm = TRUE
        )
      }
    )
    ci <- apply(
      X = draws,
      MARGIN = 2,
      FUN = function(x) {
        if (all(is.na(x))) {
          rep(
            x = NA_real_,
            times = length(probs)
          )
        } else {
          stats::quantile(
            x = x,
            probs = probs,
            na.rm = TRUE,
            names = FALSE
          )
        }
      }
    )
    ci <- t(ci)
    out <- cbind(
      est = est,
      se = se,
      R = r,
      ci
    )
    colnames(out) <- c(
      "est",
      "se",
      "R",
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
    "beta_vec",
    "kappa_vec",
    "phi_vec",
    "omega_vec",
    "psi_vec",
    "mu_x_vec",
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
  if (length(target_names) == 0L) {
    stop(
      paste(
        "No Monte Carlo confidence interval targets were found.",
        "Check that the fitted OpenMx model contains the expected",
        "matrices or algebras."
      ),
      call. = FALSE
    )
  }
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
  if (length(missing_pars) > 0L) {
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
  names(draws_list) <- target_names
  if (is.null(ncores)) {
    ncores <- 1L
  } else {
    ncores <- as.integer(ncores)
    if (is.na(ncores) || ncores < 1L) {
      ncores <- 1L
    }
    ncores <- min(
      ncores,
      parallel::detectCores(),
      nrep
    )
  }
  if (ncores > 1L) {
    # nocov start
    threads <- OpenMx::mxOption(
      key = "Number of Threads"
    )
    on.exit(
      OpenMx::mxOption(
        key = "Number of Threads",
        value = threads
      ),
      add = TRUE
    )
    OpenMx::mxOption(
      key = "Number of Threads",
      value = 1L
    )
    os_type <- Sys.info()["sysname"]
    if (os_type == "Darwin") {
      fork <- TRUE
    } else if (os_type == "Linux") {
      fork <- TRUE
    } else {
      fork <- FALSE
    }
    # nocov end
  } else {
    fork <- FALSE
  }
  .OneDraw <- function(i) {
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
      return(
        list(
          ok = FALSE,
          values = NULL
        )
      )
    }
    values <- lapply(
      X = target_names,
      FUN = function(name) {
        tryCatch(
          expr = .EvalByNameVec(
            name = name,
            model = model_i
          ),
          error = function(e) {
            rep(
              x = NA_real_,
              times = length(est_list[[name]])
            )
          }
        )
      }
    )
    names(values) <- target_names
    list(
      ok = TRUE,
      values = values
    )
  }
  if (ncores > 1L) {
    # nocov start
    if (fork) {
      draw_results <- parallel::mclapply(
        X = seq_len(nrep),
        FUN = .OneDraw,
        mc.cores = ncores
      )
    } else {
      cl <- parallel::makeCluster(ncores)
      on.exit(
        parallel::stopCluster(cl = cl),
        add = TRUE
      )
      draw_results <- parallel::parLapply(
        cl = cl,
        X = seq_len(nrep),
        fun = .OneDraw
      )
    }
    # nocov end
  } else {
    draw_results <- lapply(
      X = seq_len(nrep),
      FUN = .OneDraw
    )
  }
  ok <- vapply(
    X = draw_results,
    FUN = function(x) {
      isTRUE(x$ok)
    },
    FUN.VALUE = logical(1)
  )
  if (!any(ok)) {
    stop(
      "No valid Monte Carlo draws were generated.",
      call. = FALSE
    )
  }
  for (j in seq_along(target_names)) {
    target_name <- target_names[j]

    for (i in seq_len(nrep)) {
      if (ok[i]) {
        draws_list[[j]][i, ] <- draw_results[[i]]$values[[target_name]]
      }
    }
  }
  draws_list <- lapply(
    X = draws_list,
    FUN = function(x) {
      x[
        ok, ,
        drop = FALSE
      ]
    }
  )
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
  if ("mu_x_vec" %in% target_names) {
    out$x0 <- .CIMC(
      est = est_list$mu_x_vec,
      draws = draws_list$mu_x_vec,
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
  if (length(ie_xyz) > 0L) {
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
  if (length(out) == 0L) {
    stop(
      "No Monte Carlo confidence interval results were generated.",
      call. = FALSE
    )
  }
  out <- do.call(
    what = "rbind",
    args = out
  )
  storage.mode(out) <- "numeric"
  out
}
