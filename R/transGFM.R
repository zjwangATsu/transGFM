# ============================================================================
# transGFM Package - Complete Function Library
#
# Transfer Learning for Generalized Factor Models
#
# Main Functions:
# 1. transGFM - Single source transfer learning
# 2. transGFM_multi - Multiple source transfer learning
# 3. source_potential - Identify potential sources based on rank comparison
# 4. source_detection - Detect positive/negative transfer sources using ratio method
# 5. identify - Factor decomposition
# 6. ic_criterion - Information criterion (IC1/IC2) for rank selection
# ============================================================================
#' @importFrom stats sd rnorm runif rpois
NULL

# ============================================================================
# Part 1: Likelihood Functions for Different Data Types
# ============================================================================

lik <- function(M, data, Omega, data_type = "continuous", sigma2 = 1) {
  if (data_type == "continuous") {
    return(-sum(Omega * (data - M)^2) / (2 * sigma2))
  } else if (data_type == "count") {
    # Poisson likelihood
    M_clipped <- pmin(M, 500)  # Avoid numerical overflow
    lambda <- exp(M_clipped)
    return(sum(Omega * (data * M_clipped - lambda)))
  } else if (data_type == "binary") {
    # Binary/Binomial likelihood
    prob <- 1 / (1 + exp(-M))
    return(sum(Omega * data * M) - sum(Omega * log(1 + exp(M))))
  } else {
    stop("Unknown data_type. Use 'continuous', 'count', or 'binary'")
  }
}

lik.row <- function(M, data, Omega, data_type = "continuous", sigma2 = 1) {
  if (data_type == "continuous") {
    return(-rowSums(Omega * (data - M)^2) / (2 * sigma2))
  } else if (data_type == "count") {
    M_clipped <- pmin(M, 500)
    lambda <- exp(M_clipped)
    return(rowSums(Omega * (data * M_clipped - lambda)))
  } else if (data_type == "binary") {
    return(rowSums(Omega * data * M) - rowSums(Omega * log(1 + exp(M))))
  }
}

lik.col <- function(M, data, Omega, data_type = "continuous", sigma2 = 1) {
  if (data_type == "continuous") {
    return(-colSums(Omega * (data - M)^2) / (2 * sigma2))
  } else if (data_type == "count") {
    M_clipped <- pmin(M, 500)
    lambda <- exp(M_clipped)
    return(colSums(Omega * (data * M_clipped - lambda)))
  } else if (data_type == "binary") {
    return(colSums(Omega * data * M) - colSums(Omega * log(1 + exp(M))))
  }
}

# ============================================================================
# Part 2: Gradient Functions for Different Data Types
# ============================================================================

grad.Theta <- function(data, Omega, A0, Theta0, data_type = "continuous") {
  M0 <- Theta0 %*% t(A0)

  if (data_type == "continuous") {
    grad.M <- (data - M0) * Omega
  } else if (data_type == "count") {
    M0_clipped <- pmin(M0, 500)
    lambda0 <- exp(M0_clipped)
    grad.M <- (data - lambda0) * Omega
  } else if (data_type == "binary") {
    prob0 <- 1 / (1 + exp(-M0))
    grad.M <- (data - prob0) * Omega
  }

  grad.Theta <- grad.M %*% A0
  return(grad.Theta)
}

grad2.Theta <- function(data, Omega, A0, Theta0, data_type = "continuous") {
  M0 <- Theta0 %*% t(A0)

  if (data_type == "continuous") {
    grad2.M <- Omega
  } else if (data_type == "count") {
    M0_clipped <- pmin(M0, 500)
    lambda0 <- exp(M0_clipped)
    grad2.M <- lambda0 * Omega
  } else if (data_type == "binary") {
    prob0 <- 1 / (1 + exp(-M0))
    grad2.M <- prob0 * (1 - prob0) * Omega
  }

  grad2.Theta <- grad2.M %*% (A0^2)
  grad2.Theta[grad2.Theta < 1e-10] <- 1e-10
  return(grad2.Theta)
}

grad.A <- function(data, Omega, A0, Theta0, data_type = "continuous") {
  M0 <- Theta0 %*% t(A0)

  if (data_type == "continuous") {
    grad.M <- (data - M0) * Omega
  } else if (data_type == "count") {
    M0_clipped <- pmin(M0, 500)
    lambda0 <- exp(M0_clipped)
    grad.M <- (data - lambda0) * Omega
  } else if (data_type == "binary") {
    prob0 <- 1 / (1 + exp(-M0))
    grad.M <- (data - prob0) * Omega
  }

  grad.A <- t(grad.M) %*% Theta0
  return(grad.A)
}

grad2.A <- function(data, Omega, A0, Theta0, data_type = "continuous") {
  M0 <- Theta0 %*% t(A0)

  if (data_type == "continuous") {
    grad2.M <- Omega
  } else if (data_type == "count") {
    M0_clipped <- pmin(M0, 500)
    lambda0 <- exp(M0_clipped)
    grad2.M <- lambda0 * Omega
  } else if (data_type == "binary") {
    prob0 <- 1 / (1 + exp(-M0))
    grad2.M <- prob0 * (1 - prob0) * Omega
  }

  grad2.A <- t(grad2.M) %*% (Theta0^2)
  grad2.A[grad2.A < 1e-10] <- 1e-10
  return(grad2.A)
}

# ============================================================================
# Part 3: Projection and Search Functions
# ============================================================================

proj <- function(Theta1, C) {
  Theta1.norm <- sqrt(rowSums(Theta1^2))
  exceed_idx <- Theta1.norm > C
  if (any(exceed_idx)) {
    Theta1[exceed_idx, ] <- Theta1[exceed_idx, ] / Theta1.norm[exceed_idx] * C
  }
  return(Theta1)
}

search.Theta <- function(grad, grad2, data, Omega, A0, Theta0,
                         data_type, step, times, sigma2 = 1) {
  n <- nrow(Theta0)
  lik.row0 <- lik.row(Theta0 %*% t(A0), data, Omega, data_type, sigma2)
  step.vec <- step * rep(1, n)
  Theta1 <- Theta0 + step.vec * grad / grad2
  lik.row1 <- lik.row(Theta1 %*% t(A0), data, Omega, data_type, sigma2)

  z <- 0
  while (min(lik.row1 - lik.row0) < 0 & z < times) {
    step.vec[lik.row1 - lik.row0 < 0] <- step.vec[lik.row1 - lik.row0 < 0] / 2
    z <- z + 1
    Theta1 <- Theta0 + step.vec * grad / grad2
    lik.row1 <- lik.row(Theta1 %*% t(A0), data, Omega, data_type, sigma2)
  }
  Theta1[lik.row1 - lik.row0 < 0, ] <- Theta0[lik.row1 - lik.row0 < 0, ]
  return(Theta1)
}

search.A <- function(grad, grad2, data, Omega, A0, Theta0,
                     data_type, step, times, sigma2 = 1) {
  p <- nrow(A0)
  lik.col0 <- lik.col(Theta0 %*% t(A0), data, Omega, data_type, sigma2)
  step.vec <- step * rep(1, p)
  A1 <- A0 + step.vec * grad / grad2
  lik.col1 <- lik.col(Theta0 %*% t(A1), data, Omega, data_type, sigma2)

  z <- 0
  while (min(lik.col1 - lik.col0) < 0 & z < times) {
    step.vec[lik.col1 - lik.col0 < 0] <- step.vec[lik.col1 - lik.col0 < 0] / 2
    z <- z + 1
    A1 <- A0 + step.vec * grad / grad2
    lik.col1 <- lik.col(Theta0 %*% t(A1), data, Omega, data_type, sigma2)
  }
  A1[lik.col1 - lik.col0 < 0, ] <- A0[lik.col1 - lik.col0 < 0, ]
  return(A1)
}

search.Theta.ball <- function(grad, grad2, data, Omega, A0, Theta0,
                              data_type, C, step, times, sigma2 = 1) {
  n <- nrow(Theta0)
  lik.row0 <- lik.row(Theta0 %*% t(A0), data, Omega, data_type, sigma2)
  step.vec <- step * rep(1, n)
  Theta1 <- proj(Theta0 + step.vec * grad / grad2, C)
  lik.row1 <- lik.row(Theta1 %*% t(A0), data, Omega, data_type, sigma2)

  z <- 0
  while (min(lik.row1 - lik.row0) < 0 & z < times) {
    step.vec[lik.row1 - lik.row0 < 0] <- step.vec[lik.row1 - lik.row0 < 0] / 2
    z <- z + 1
    Theta1 <- proj(Theta0 + step.vec * grad / grad2, C)
    lik.row1 <- lik.row(Theta1 %*% t(A0), data, Omega, data_type, sigma2)
  }
  Theta1[lik.row1 - lik.row0 < 0, ] <- Theta0[lik.row1 - lik.row0 < 0, ]
  return(Theta1)
}

search.A.ball <- function(grad, grad2, data, Omega, A0, Theta0,
                          data_type, C, step, times, sigma2 = 1) {
  p <- nrow(A0)
  lik.col0 <- lik.col(Theta0 %*% t(A0), data, Omega, data_type, sigma2)
  step.vec <- step * rep(1, p)
  A1 <- proj(A0 + step.vec * grad / grad2, C)
  lik.col1 <- lik.col(Theta0 %*% t(A1), data, Omega, data_type, sigma2)

  z <- 0
  while (min(lik.col1 - lik.col0) < 0 & z < times) {
    step.vec[lik.col1 - lik.col0 < 0] <- step.vec[lik.col1 - lik.col0 < 0] / 2
    z <- z + 1
    A1 <- proj(A0 + step.vec * grad / grad2, C)
    lik.col1 <- lik.col(Theta0 %*% t(A1), data, Omega, data_type, sigma2)
  }
  A1[lik.col1 - lik.col0 < 0, ] <- A0[lik.col1 - lik.col0 < 0, ]
  return(A1)
}

# ============================================================================
# Part 4: Core Algorithm - CJMLE
# ============================================================================

CJMLE <- function(Theta0, A0, r, data, Omega, C,
                  data_type = "continuous", method = "newton",
                  sigma2 = 1, tol = 1e-4, max_iter = 30, verbose = FALSE) {

  # Initialize
  Theta <- Theta0
  A <- A0
  n <- nrow(data)
  p <- ncol(data)

  # Store history
  obj_hist <- numeric(max_iter)

  # Main iteration loop
  for (iter in 1:max_iter) {
    Theta_old <- Theta
    A_old <- A

    # Update Theta
    grad_theta <- grad.Theta(data, Omega, A, Theta, data_type)

    if (method == "newton") {
      grad2_theta <- grad2.Theta(data, Omega, A, Theta, data_type)
      Theta <- search.Theta.ball(grad_theta, grad2_theta, data, Omega, A, Theta,
                                 data_type, C, step = 1, times = 10, sigma2)
    } else {
      Theta <- proj(Theta + 0.01 * grad_theta, C)
    }

    # Update A
    grad_a <- grad.A(data, Omega, A, Theta, data_type)

    if (method == "newton") {
      grad2_a <- grad2.A(data, Omega, A, Theta, data_type)
      A <- search.A.ball(grad_a, grad2_a, data, Omega, A, Theta,
                         data_type, C, step = 1, times = 10, sigma2)
    } else {
      A <- proj(A + 0.01 * grad_a, C)
    }

    # Compute objective
    M <- Theta %*% t(A)
    obj_hist[iter] <- lik(M, data, Omega, data_type, sigma2)

    # Check convergence
    theta_change <- sqrt(sum((Theta - Theta_old)^2)) / (sqrt(sum(Theta_old^2)) + 1e-10)
    a_change <- sqrt(sum((A - A_old)^2)) / (sqrt(sum(A_old^2)) + 1e-10)

    if (verbose) {
      cat(sprintf("Iter %d: Obj = %.4f, Theta change = %.6f, A change = %.6f\n",
                  iter, obj_hist[iter], theta_change, a_change))
    }

    if (max(theta_change, a_change) < tol && iter > 5) {
      if (verbose) cat(sprintf("Converged at iteration %d\n", iter))
      obj_hist <- obj_hist[1:iter]
      break
    }
  }

  return(list(
    M = Theta %*% t(A),
    Theta = Theta,
    A = A,
    iterations = iter,
    objective = obj_hist[length(obj_hist)],
    objective_history = obj_hist,
    converged = (iter < max_iter)
  ))
}

# ============================================================================
# Part 5: Refinement Algorithm
# ============================================================================

refi.nosp <- function(M0, r, data, Omega, C, data_type = "continuous",
                      sigma2 = 1, tol = 1e-4, max_iter = 30, verbose = FALSE) {

  # Initial decomposition
  svd_M0 <- svd(M0)
  d_sqrt <- sqrt(abs(svd_M0$d[1:r]))
  Theta <- svd_M0$u[, 1:r, drop = FALSE] %*% diag(d_sqrt, nrow = r)
  A <- svd_M0$v[, 1:r, drop = FALSE] %*% diag(d_sqrt, nrow = r)

  n <- nrow(data)
  p <- ncol(data)
  obj_hist <- numeric(max_iter)

  # Main iteration loop
  for (iter in 1:max_iter) {
    Theta_old <- Theta
    A_old <- A

    # Update Theta
    grad_theta <- grad.Theta(data, Omega, A, Theta, data_type)
    grad2_theta <- grad2.Theta(data, Omega, A, Theta, data_type)
    Theta <- search.Theta(grad_theta, grad2_theta, data, Omega, A, Theta,
                          data_type, step = 1, times = 10, sigma2)

    # Update A
    grad_a <- grad.A(data, Omega, A, Theta, data_type)
    grad2_a <- grad2.A(data, Omega, A, Theta, data_type)
    A <- search.A(grad_a, grad2_a, data, Omega, A, Theta,
                  data_type, step = 1, times = 10, sigma2)

    # Compute objective
    M <- Theta %*% t(A)
    obj_hist[iter] <- lik(M, data, Omega, data_type, sigma2)

    # Check convergence
    theta_change <- sqrt(sum((Theta - Theta_old)^2)) / (sqrt(sum(Theta_old^2)) + 1e-10)
    a_change <- sqrt(sum((A - A_old)^2)) / (sqrt(sum(A_old^2)) + 1e-10)

    if (verbose) {
      cat(sprintf("Iter %d: Obj = %.4f, Theta change = %.6f, A change = %.6f\n",
                  iter, obj_hist[iter], theta_change, a_change))
    }

    if (max(theta_change, a_change) < tol && iter > 5) {
      if (verbose) cat(sprintf("Converged at iteration %d\n", iter))
      obj_hist <- obj_hist[1:iter]
      break
    }
  }

  return(list(
    M = Theta %*% t(A),
    Theta = Theta,
    A = A,
    iterations = iter,
    objective = obj_hist[length(obj_hist)],
    objective_history = obj_hist,
    converged = (iter < max_iter)
  ))
}

# ============================================================================
# Part 6: Nuclear Penalized MLE Helper Functions
# ============================================================================

nuclear_norm <- function(M) {
  s <- svd(M, nu = 0, nv = 0)
  return(sum(s$d))
}

element_wise_projection <- function(M, Cd) {
  if (!is.finite(Cd)) return(M)
  M[M > Cd] <- Cd
  M[M < -Cd] <- -Cd
  return(M)
}

soft_threshold_svd_with_projection <- function(Z, lambda, Cd) {
  svd_Z <- svd(Z)
  d_thresh <- pmax(svd_Z$d - lambda, 0)

  if (sum(d_thresh > 0) == 0) {
    M_new <- matrix(0, nrow(Z), ncol(Z))
  } else {
    r_eff <- sum(d_thresh > 0)
    M_new <- svd_Z$u[, 1:r_eff, drop = FALSE] %*%
      diag(d_thresh[1:r_eff], nrow = r_eff) %*%
      t(svd_Z$v[, 1:r_eff, drop = FALSE])
  }

  if (is.finite(Cd)) {
    M_new <- element_wise_projection(M_new, Cd)
  }

  return(M_new)
}

compute_loglikelihood <- function(X, C, M, Delta, family, sigma2 = 1) {
  M_total <- M + Delta

  if (family == "continuous") {
    return(-sum(C * (X - M_total)^2) / (2 * sigma2))
  } else if (family == "count") {
    M_total_clipped <- pmin(M_total, 500)
    lambda <- exp(M_total_clipped)
    return(sum(C * (X * M_total_clipped - lambda)))
  } else if (family == "binary") {
    return(sum(C * X * M_total) - sum(C * log(1 + exp(M_total))))
  } else {
    stop("Unknown family type")
  }
}

compute_gradient_improved <- function(X, C, M, Delta, family, sigma2 = 1) {
  M_total <- M + Delta

  if (family == "continuous") {
    grad <- (X - M_total) * C / sigma2
  } else if (family == "count") {
    M_total_clipped <- pmin(M_total, 500)
    lambda <- exp(M_total_clipped)
    grad <- (X - lambda) * C
  } else if (family == "binary") {
    prob <- 1 / (1 + exp(-M_total))
    grad <- (X - prob) * C
  } else {
    stop("Unknown family type")
  }

  return(grad)
}

adaptive_newton_update <- function(X, C, M, Delta, family, lambda, Cd, sigma2 = 1) {
  M_total <- M + Delta
  n <- nrow(X)
  p <- ncol(X)

  # Compute gradient
  grad <- compute_gradient_improved(X, C, M, Delta, family, sigma2)

  # Compute Hessian diagonal approximation
  if (family == "continuous") {
    H_diag <- C / sigma2
  } else if (family == "count") {
    M_total_clipped <- pmin(M_total, 500)
    lambda_val <- exp(M_total_clipped)
    H_diag <- lambda_val * C
  } else if (family == "binary") {
    prob <- 1 / (1 + exp(-M_total))
    H_diag <- prob * (1 - prob) * C
  }

  H_diag <- pmax(H_diag, 1e-10)

  # Newton direction
  newton_dir <- grad / H_diag

  # Line search
  step <- 1.0
  max_backtrack <- 20

  current_obj <- -compute_loglikelihood(X, C, M, Delta, family, sigma2) +
    lambda * nuclear_norm(Delta)

  for (bt in 1:max_backtrack) {
    Z <- Delta + step * newton_dir
    Delta_new <- soft_threshold_svd_with_projection(Z, step * lambda, Cd)

    new_obj <- -compute_loglikelihood(X, C, M, Delta_new, family, sigma2) +
      lambda * nuclear_norm(Delta_new)

    if (is.finite(new_obj) && new_obj < current_obj) {
      break
    }

    step <- step * 0.5
  }

  return(list(Delta_new = Delta_new, step = step))
}

# ============================================================================
# Part 7: Nuclear Penalized MLE
# ============================================================================

nuclear_penalized_mle_improved <- function(X, C, M, family = "count",
                                           lambda = 1, Cd = Inf, sigma2 = 1,
                                           max_iter = 100, tol = 1e-5,
                                           verbose = FALSE, warm_start = TRUE,
                                           use_newton = TRUE) {

  n <- nrow(X)
  p <- ncol(X)

  # Initialize Delta
  if (warm_start) {
    # Warm start using CJMLE-based initialization
    residual <- X - M
    svd_residual <- svd(residual)
    rank_init <- min(3, sum(svd_residual$d > 1e-2))

    if (rank_init > 0) {
      Delta <- svd_residual$u[, 1:rank_init, drop = FALSE] %*%
        diag(svd_residual$d[1:rank_init], nrow = rank_init) %*%
        t(svd_residual$v[, 1:rank_init, drop = FALSE])
    } else {
      Delta <- residual * 0.01
    }
  } else {
    Delta <- matrix(0, n, p)
  }

  # Apply initial projection
  if (is.finite(Cd)) {
    Delta <- element_wise_projection(Delta, Cd)
  }

  # Store history
  obj_values <- numeric(max_iter)
  rank_values <- numeric(max_iter)

  # Main optimization loop
  for (iter in 1:max_iter) {
    Delta_old <- Delta

    if (use_newton) {
      # Newton update
      update_result <- adaptive_newton_update(X, C, M, Delta, family, lambda, Cd, sigma2)
      Delta <- update_result$Delta_new
    } else {
      # Gradient descent
      grad <- compute_gradient_improved(X, C, M, Delta, family, sigma2)
      step_size <- 0.1
      Z <- Delta + step_size * grad
      Delta <- soft_threshold_svd_with_projection(Z, step_size * lambda, Cd)
    }

    # Check for invalid values
    if (any(is.na(Delta) | is.infinite(Delta))) {
      warning("Delta contains invalid values, reverting")
      Delta <- Delta_old
      break
    }

    # Compute objective
    current_loglik <- compute_loglikelihood(X, C, M, Delta, family, sigma2)
    current_penalty <- lambda * nuclear_norm(Delta)
    obj_values[iter] <- -current_loglik + current_penalty

    # Record rank
    svd_delta <- svd(Delta, nu = 0, nv = 0)
    rank_values[iter] <- sum(svd_delta$d > 1e-3)

    # Check convergence
    change <- sqrt(sum((Delta - Delta_old)^2)) / (sqrt(sum(Delta_old^2)) + 1e-10)

    if (verbose && (iter <= 10 || iter %% 50 == 0)) {
      cat(sprintf("Iter %d: Obj = %.4f, Rank = %d, Change = %.6f\n",
                  iter, obj_values[iter], rank_values[iter], change))
    }

    if (change < tol && iter > 10) {
      if (verbose) cat(sprintf("Converged at iteration %d\n", iter))
      obj_values <- obj_values[1:iter]
      rank_values <- rank_values[1:iter]
      break
    }
  }

  # Final projection
  Delta <- element_wise_projection(Delta, Cd)

  # Final diagnostics
  final_svd <- svd(Delta, nu = 0, nv = 0)
  final_rank <- sum(final_svd$d > 1e-3)
  final_loglik <- compute_loglikelihood(X, C, M, Delta, family, sigma2)
  final_penalty <- lambda * nuclear_norm(Delta)

  return(list(
    Delta = Delta,
    loglikelihood = final_loglik,
    penalized_loglikelihood = final_loglik - final_penalty,
    iterations = iter,
    nuclear_norm = nuclear_norm(Delta),
    rank = final_rank,
    singular_values = final_svd$d,
    objective_values = obj_values,
    rank_history = rank_values,
    converged = (iter < max_iter),
    lambda = lambda,
    sigma2 = sigma2
  ))
}

# ============================================================================
# Part 8: Cross-Validation for Lambda Selection
# ============================================================================

cv_nuclear_penalized_mle <- function(X, C, M, family = "count",
                                     lambda_seq = seq(0, 10, by = 1),
                                     K = 3, max_iter = 100, Cd = Inf,
                                     sigma2 = 1, verbose = FALSE) {

  n <- nrow(X)
  p <- ncol(X)
  n_lambda <- length(lambda_seq)

  # Create CV folds
  fold_ids <- sample(rep(1:K, length.out = sum(C)))
  cv_errors <- matrix(0, nrow = n_lambda, ncol = K)

  for (k in 1:K) {
    if (verbose) cat(sprintf("CV Fold %d/%d\n", k, K))

    # Create train/test split
    C_train <- C
    C_test <- matrix(0, n, p)

    obs_idx <- which(C == 1)
    test_obs <- obs_idx[fold_ids == k]
    test_rows <- ((test_obs - 1) %% n) + 1
    test_cols <- ((test_obs - 1) %/% n) + 1

    for (i in seq_along(test_obs)) {
      C_train[test_rows[i], test_cols[i]] <- 0
      C_test[test_rows[i], test_cols[i]] <- 1
    }

    # Try each lambda
    for (lambda_idx in 1:n_lambda) {
      lambda <- lambda_seq[lambda_idx]

      # Fit model on training set
      result <- nuclear_penalized_mle_improved(
        X = X, C = C_train, M = M, family = family,
        lambda = lambda, Cd = Cd, sigma2 = sigma2,
        max_iter = max_iter, verbose = FALSE,
        warm_start = TRUE, use_newton = TRUE
      )

      # Compute test error
      M_pred <- M + result$Delta
      test_loglik <- compute_loglikelihood(X, C_test, matrix(0, n, p), M_pred, family, sigma2)
      cv_errors[lambda_idx, k] <- -test_loglik / sum(C_test)
    }
  }

  # Average CV errors
  mean_cv_errors <- rowMeans(cv_errors)
  se_cv_errors <- apply(cv_errors, 1, sd) / sqrt(K)

  # Select optimal lambda (minimum CV error)
  opt_idx <- which.min(mean_cv_errors)
  opt_lambda <- lambda_seq[opt_idx]

  if (verbose) {
    cat(sprintf("\nOptimal lambda = %.2f (CV error = %.4f)\n",
                opt_lambda, mean_cv_errors[opt_idx]))
  }

  return(list(
    optimal_lambda = opt_lambda,
    lambda_seq = lambda_seq,
    mean_cv_errors = mean_cv_errors,
    se_cv_errors = se_cv_errors,
    cv_errors_all = cv_errors
  ))
}

# ============================================================================
# Part 9: Utility Functions
# ============================================================================

#' Identify factor decomposition via SVD
#'
#' @param M Matrix to decompose
#' @param r Number of factors
#' @return List with F (row factors) and B (column factors)
#'
#' @examples
#' # Generate Poisson data
#' set.seed(123)
#' n0 <- 50; p0 <- 50; r <- 2
#' F_true <- matrix(runif(n0 * r, min = -2, max = 2), n0, r)
#' B_true <- matrix(runif(p0 * r, min = -2, max = 2), p0, r)
#' F_true <- F_true / sqrt(r)
#' B_true <- B_true / sqrt(r)
#' M_true <- F_true %*% t(B_true)
#'
#' # Decompose using identify
#' result <- identify(M_true, r = 2)
#' F_hat <- result$F
#' B_hat <- result$B
#'
#' # Check reconstruction
#' M_reconstructed <- F_hat %*% t(B_hat)
#' print(max(abs(M_reconstructed - M_true)))  # Should be very small
#'
#' @export
identify <- function(M, r) {
  svd_M <- svd(M)
  d_sqrt <- sqrt(abs(svd_M$d[1:r]))
  F_mat <- svd_M$u[, 1:r, drop = FALSE] %*% diag(d_sqrt, nrow = r)
  B_mat <- svd_M$v[, 1:r, drop = FALSE] %*% diag(d_sqrt, nrow = r)

  return(list(F = F_mat, B = B_mat))
}

best_rank_r <- function(M, r) {
  svd_M <- svd(M)
  if (r == 0) {
    return(matrix(0, nrow(M), ncol(M)))
  }

  M_r <- svd_M$u[, 1:r, drop = FALSE] %*%
    diag(svd_M$d[1:r], nrow = r) %*%
    t(svd_M$v[, 1:r, drop = FALSE])

  return(M_r)
}
#' Calculate relative error between estimated and true matrices
#'
#' @param M_hat Estimated matrix
#' @param M_true True matrix
#' @return Relative Frobenius norm error
#' @export
#'
#' @examples
#' M_true <- matrix(1:9, 3, 3)
#' M_hat <- M_true + matrix(rnorm(9, 0, 0.1), 3, 3)
#' relative_error(M_hat, M_true)
relative_error <- function(M_hat, M_true) {
  return(sqrt(sum((M_hat - M_true)^2)) / sqrt(sum(M_true^2)))
}

calculate_projection_constants <- function(r, p, upper=2) {
  # CJMLE stage: C = ceiling(sqrt(r * 4))
  C_cjmle <- ceiling(sqrt(r * upper* upper))

  # Refinement stage: C2 = 2 * sqrt(r/p)
  C_refine <- 2 * sqrt(r/p)

  return(list(C = C_cjmle, C2 = C_refine))
}

# ============================================================================
# Part 10: Main Transfer Learning Functions
# ============================================================================

#' Single source transfer learning for generalized factor models
#'
#' @param source_data Source data matrix (may contain missing values coded as NA)
#' @param target_data Target data matrix (complete)
#' @param r Number of factors
#' @param data_type Type of data: "continuous", "count", or "binary"
#' @param lambda_seq Sequence of lambda values for CV (default: seq(0, 10, by = 1))
#' @param K_cv Number of CV folds (default: 3)
#' @param sigma2 Variance parameter for continuous data (default: 1)
#' @param max_iter_cjmle Maximum iterations for CJMLE (default: 30)
#' @param max_iter_refine Maximum iterations for refinement (default: 30)
#' @param max_iter_nuclear Maximum iterations for nuclear MLE (default: 100)
#' @param verbose Print progress information (default: FALSE)
#' @return List containing final estimate M_trans and intermediate results
#'
#' @examples
#' # Generate Poisson data
#' set.seed(2025)
#'
#' # Source data (100 x 100 with 10% missing)
#' n1 <- 100; p1 <- 100; r <- 2
#' F_source <- matrix(runif(n1 * r, min = -2, max = 2), n1, r)
#' B_source <- matrix(runif(p1 * r, min = -2, max = 2), p1, r)
#' M_source <- F_source %*% t(B_source)
#' lambda_source <- exp(M_source)
#' X_source <- matrix(rpois(n1 * p1, as.vector(lambda_source)), n1, p1)
#'
#' # Add 10% missing values to source
#' n_missing <- floor(n1 * p1 * 0.1)
#' missing_idx <- sample(n1 * p1, n_missing)
#' X_source[missing_idx] <- NA
#'
#' # Target data (50 x 50, complete)
#' n0 <- 50; p0 <- 50
#' M_target_true <- M_source[1:n0, 1:p0]
#' lambda_target <- exp(M_target_true)
#' X_target <- matrix(rpois(n0 * p0, as.vector(lambda_target)), n0, p0)
#'
#' # Run transGFM
#' result <- transGFM(
#'   source_data = X_source,
#'   target_data = X_target,
#'   r = 2,
#'   data_type = "count",
#'   lambda_seq = seq(0, 5, by = 1),
#'   K_cv = 3,
#'   verbose = FALSE
#' )
#'
#' # Check results
#' print(paste("Optimal lambda:", result$optimal_lambda))
#' print(paste("Final relative error:",
#'             relative_error(result$M_trans, M_target_true)))
#'
#' @export
transGFM <- function(source_data, target_data, r, data_type = "count",
                     lambda_seq = seq(0, 10, by = 1), K_cv = 3,
                     sigma2 = 1, max_iter_cjmle = 30, max_iter_refine = 30,
                     max_iter_nuclear = 30, verbose = FALSE) {

  if (verbose) cat("=== Starting transGFM (Single Source) ===\n\n")

  # Get dimensions
  n1 <- nrow(source_data)
  p1 <- ncol(source_data)
  n0 <- nrow(target_data)
  p0 <- ncol(target_data)

  # Check dimension compatibility
  if (n0 > n1 || p0 > p1) {
    stop("Target dimensions must be smaller than or equal to source dimensions")
  }

  # Create observation indicator matrix for source
  Omega_source <- matrix(1, nrow = n1, ncol = p1)
  Omega_source[is.na(source_data)] <- 0
  source_data_filled <- source_data
  source_data_filled[is.na(source_data_filled)] <- 0

  # Calculate projection constants
  proj_consts <- calculate_projection_constants(r, p1, upper = 2)
  C <- proj_consts$C
  C2 <- proj_consts$C2
  C3 <- 2 * C^2

  if (verbose) {
    cat(sprintf("Source dimensions: %d x %d\n", n1, p1))
    cat(sprintf("Target dimensions: %d x %d\n", n0, p0))
    cat(sprintf("Missing rate in source: %.2f%%\n",
                (1 - mean(Omega_source)) * 100))
    cat(sprintf("Projection constants: C=%.2f, C2=%.2f, C3=%.2f\n\n", C, C2, C3))
  }

  # Step 1: Initialize factors
  if (verbose) cat("Step 1: Initializing factors...\n")
  if (data_type == "count") {
    log_data <- log(source_data_filled + 1)
  } else if (data_type == "binary") {
    log_data <- source_data_filled
  } else {
    log_data <- source_data_filled
  }

  svd_init <- svd(log_data)
  d_sqrt <- sqrt(abs(svd_init$d[1:r]))
  F_init <- svd_init$u[, 1:r, drop = FALSE] %*% diag(d_sqrt, nrow = r)
  B_init <- svd_init$v[, 1:r, drop = FALSE] %*% diag(d_sqrt, nrow = r)

  # Step 2: CJMLE on source data
  if (verbose) cat("Step 2: Running CJMLE on source data...\n")
  result_cjmle <- CJMLE(
    Theta0 = F_init, A0 = B_init, r = r,
    data = source_data_filled, Omega = Omega_source, C = C,
    data_type = data_type, method = "newton",
    sigma2 = sigma2, max_iter = max_iter_cjmle, verbose = verbose
  )
  M1_hat <- result_cjmle$M

  if (verbose) {
    cat(sprintf("CJMLE completed in %d iterations\n\n", result_cjmle$iterations))
  }

  # Step 3: Refinement
  if (verbose) cat("Step 3: Refining estimate...\n")
  result_refine <- refi.nosp(
    M0 = M1_hat, r = r, data = source_data_filled,
    Omega = Omega_source, C = C2, data_type = data_type,
    sigma2 = sigma2, max_iter = max_iter_refine, verbose = verbose
  )
  M1_refined <- result_refine$M

  if (verbose) {
    cat(sprintf("Refinement completed in %d iterations\n\n",
                result_refine$iterations))
  }

  # Step 4: Extract target region
  if (verbose) cat("Step 4: Extracting target region...\n")
  M10_refined <- M1_refined[1:n0, 1:p0]

  # Step 5: Cross-validation for lambda
  if (verbose) cat("Step 5: Cross-validation for lambda selection...\n")
  cv_result <- cv_nuclear_penalized_mle(
    X = target_data,
    C = matrix(1, nrow = n0, ncol = p0),
    M = M10_refined,
    family = data_type,
    lambda_seq = lambda_seq,
    K = K_cv,
    max_iter = max_iter_nuclear,
    Cd = C3,
    sigma2 = sigma2,
    verbose = verbose
  )
  opt_lambda <- cv_result$optimal_lambda

  if (verbose) {
    cat(sprintf("Optimal lambda selected: %.2f\n\n", opt_lambda))
  }

  # Step 6: Debiased estimation with nuclear penalty
  if (verbose) cat("Step 6: Debiased estimation with nuclear penalty...\n")
  result_nuclear <- nuclear_penalized_mle_improved(
    X = target_data,
    C = matrix(1, nrow = n0, ncol = p0),
    M = M10_refined,
    family = data_type,
    lambda = opt_lambda,
    Cd = C3,
    sigma2 = sigma2,
    max_iter = max_iter_nuclear,
    verbose = verbose,
    warm_start = TRUE,
    use_newton = TRUE
  )

  if (verbose) {
    cat(sprintf("Nuclear MLE completed in %d iterations\n\n",
                result_nuclear$iterations))
  }

  # Step 7: Final estimate and rank-r truncation
  if (verbose) cat("Step 7: Computing final estimate...\n")
  M_debiased <- M10_refined + result_nuclear$Delta
  M_trans <- best_rank_r(M_debiased, r = r)

  if (verbose) {
    cat("=== transGFM completed ===\n")
  }

  return(list(
    M_trans = M_trans,
    M_debiased = M_debiased,
    M10_refined = M10_refined,
    M1_refined = M1_refined,
    M1_cjmle = M1_hat,
    Delta = result_nuclear$Delta,
    optimal_lambda = opt_lambda,
    cv_result = cv_result,
    cjmle_result = result_cjmle,
    refine_result = result_refine,
    nuclear_result = result_nuclear
  ))
}

#' Multiple source transfer learning for generalized factor models
#'
#' @param source_data_list List of source data matrices (may contain missing values)
#' @param target_data Target data matrix (complete)
#' @param r Number of factors
#' @param data_type Type of data: "continuous", "count", or "binary"
#' @param method Fusion method: "AD" (Average-Debias) or "DA" (Debias-Average)
#' @param lambda_seq Sequence of lambda values for CV
#' @param K_cv Number of CV folds
#' @param sigma2 Variance parameter for continuous data
#' @param max_iter_cjmle Maximum iterations for CJMLE
#' @param max_iter_refine Maximum iterations for refinement
#' @param max_iter_nuclear Maximum iterations for nuclear MLE
#' @param verbose Print progress information
#' @return List containing final estimate and intermediate results
#'
#' @examples
#' \donttest{
#' # Generate Poisson data
#' set.seed(2025)
#'
#' # Generate 3 source datasets (100 x 100 with different missing rates)
#' n1 <- 100; p1 <- 100; r <- 2
#' source_list <- list()
#' F_s <- matrix(runif(n1 * r, min = -2, max = 2), n1, r)
#' B_s <- matrix(runif(p1 * r, min = -2, max = 2), p1, r)
#' M_s <- F_s %*% t(B_s)
#' for (s in 1:3) {
#'   X_s <- matrix(rpois(n1 * p1, exp(M_s)), n1, p1)
#'
#'   # Add missing values (10%, 12%, 14% for sources 1-3)
#'   missing_rate <- 0.1 + (s - 1) * 0.02
#'   n_missing <- floor(n1 * p1 * missing_rate)
#'   missing_idx <- sample(n1 * p1, n_missing)
#'   X_s[missing_idx] <- NA
#'
#'   source_list[[s]] <- X_s
#' }
#'
#' # Target data (50 x 50, complete)
#' n0 <- 50; p0 <- 50
#' M_target_true <- M_s[1:n0, 1:p0]
#' X_target <- matrix(rpois(n0 * p0, exp(M_target_true)), n0, p0)
#'
#' # Run transGFM_multi with AD method
#' result_AD <- transGFM_multi(
#'   source_data_list = source_list,
#'   target_data = X_target,
#'   r = 2,
#'   data_type = "count",
#'   method = "AD",
#'   lambda_seq = seq(0, 5, by = 1),
#'   K_cv = 3,
#'   verbose = FALSE
#' )
#'
#' # Run transGFM_multi with DA method
#' result_DA <- transGFM_multi(
#'   source_data_list = source_list,
#'   target_data = X_target,
#'   r = 2,
#'   data_type = "count",
#'   method = "DA",
#'   verbose = FALSE
#' )
#'
#' # Compare results
#' print(paste("AD method error:", relative_error(result_AD$M_trans, M_target_true)))
#' print(paste("DA method error:", relative_error(result_DA$M_trans, M_target_true)))
#' }
#' @export
transGFM_multi <- function(source_data_list, target_data, r, data_type = "count",
                           method = "AD", lambda_seq = seq(0, 10, by = 1),
                           K_cv = 3, sigma2 = 1, max_iter_cjmle = 30,
                           max_iter_refine = 30, max_iter_nuclear = 100,
                           verbose = FALSE) {

  if (verbose) cat("=== Starting transGFM_multi (Multiple Sources) ===\n\n")

  num_sources <- length(source_data_list)
  n0 <- nrow(target_data)
  p0 <- ncol(target_data)

  if (verbose) {
    cat(sprintf("Number of sources: %d\n", num_sources))
    cat(sprintf("Target dimensions: %d x %d\n", n0, p0))
    cat(sprintf("Fusion method: %s\n\n", method))
  }

  # Calculate projection constants (using first source dimensions)
  n1 <- nrow(source_data_list[[1]])
  p1 <- ncol(source_data_list[[1]])
  proj_consts <- calculate_projection_constants(r, p1, upper = 1)
  C <- proj_consts$C
  C2 <- proj_consts$C2
  C3 <- 2 * C^2

  # Process each source
  M_refined_list <- list()
  M10_refined_list <- list()

  for (s in 1:num_sources) {
    if (verbose) {
      cat(sprintf("--- Processing source %d/%d ---\n", s, num_sources))
    }

    source_data <- source_data_list[[s]]
    n1 <- nrow(source_data)
    p1 <- ncol(source_data)

    # Create observation indicator
    Omega_source <- matrix(1, nrow = n1, ncol = p1)
    Omega_source[is.na(source_data)] <- 0
    source_data_filled <- source_data
    source_data_filled[is.na(source_data_filled)] <- 0

    # Initialize
    if (data_type == "count") {
      log_data <- log(source_data_filled + 1)
    } else {
      log_data <- source_data_filled
    }

    svd_init <- svd(log_data)
    d_sqrt <- sqrt(abs(svd_init$d[1:r]))
    F_init <- svd_init$u[, 1:r, drop = FALSE] %*% diag(d_sqrt, nrow = r)
    B_init <- svd_init$v[, 1:r, drop = FALSE] %*% diag(d_sqrt, nrow = r)

    # CJMLE
    if (verbose) cat("  Running CJMLE...\n")
    result_cjmle <- CJMLE(
      Theta0 = F_init, A0 = B_init, r = r,
      data = source_data_filled, Omega = Omega_source, C = C,
      data_type = data_type, method = "newton",
      sigma2 = sigma2, max_iter = max_iter_cjmle, verbose = FALSE
    )

    # Refinement
    if (verbose) cat("  Running refinement...\n")
    result_refine <- refi.nosp(
      M0 = result_cjmle$M, r = r, data = source_data_filled,
      Omega = Omega_source, C = C2, data_type = data_type,
      sigma2 = sigma2, max_iter = max_iter_refine, verbose = FALSE
    )

    M_refined_list[[s]] <- result_refine$M
    M10_refined_list[[s]] <- result_refine$M[1:n0, 1:p0]

    if (verbose) cat("  Source processing completed\n\n")
  }

  # Fusion based on method
  if (method == "AD") {
    # Method 1: Average then Debias
    if (verbose) cat("Using AD method: Average -> Debias\n")

    # Average refined estimates
    M10_avg <- Reduce("+", M10_refined_list) / num_sources

    # CV for lambda
    if (verbose) cat("Cross-validation for lambda...\n")
    cv_result <- cv_nuclear_penalized_mle(
      X = target_data,
      C = matrix(1, nrow = n0, ncol = p0),
      M = M10_avg,
      family = data_type,
      lambda_seq = lambda_seq,
      K = K_cv,
      max_iter = max_iter_nuclear,
      Cd = C3,
      sigma2 = sigma2,
      verbose = FALSE
    )
    opt_lambda <- cv_result$optimal_lambda

    # Debias on average
    if (verbose) cat("Debiasing averaged estimate...\n")
    result_nuclear <- nuclear_penalized_mle_improved(
      X = target_data,
      C = matrix(1, nrow = n0, ncol = p0),
      M = M10_avg,
      family = data_type,
      lambda = opt_lambda,
      Cd = C3,
      sigma2 = sigma2,
      max_iter = max_iter_nuclear,
      verbose = FALSE,
      warm_start = TRUE,
      use_newton = TRUE
    )

    M_debiased <- M10_avg + result_nuclear$Delta
    M_trans <- best_rank_r(M_debiased, r = r)

    return(list(
      M_trans = M_trans,
      M_debiased = M_debiased,
      M10_avg = M10_avg,
      M10_refined_list = M10_refined_list,
      M_refined_list = M_refined_list,
      Delta = result_nuclear$Delta,
      optimal_lambda = opt_lambda,
      cv_result = cv_result,
      nuclear_result = result_nuclear,
      method = "AD",
      num_sources = num_sources
    ))

  } else if (method == "DA") {
    # Method 2: Debias then Average
    if (verbose) cat("Using DA method: Debias -> Average\n")

    M_debiased_list <- list()

    for (s in 1:num_sources) {
      if (verbose) cat(sprintf("Debiasing source %d/%d...\n", s, num_sources))

      # CV for each source
      cv_result_s <- cv_nuclear_penalized_mle(
        X = target_data,
        C = matrix(1, nrow = n0, ncol = p0),
        M = M10_refined_list[[s]],
        family = data_type,
        lambda_seq = lambda_seq,
        K = K_cv,
        max_iter = max_iter_nuclear,
        Cd = C3,
        sigma2 = sigma2,
        verbose = FALSE
      )

      # Debias each source
      result_nuclear_s <- nuclear_penalized_mle_improved(
        X = target_data,
        C = matrix(1, nrow = n0, ncol = p0),
        M = M10_refined_list[[s]],
        family = data_type,
        lambda = cv_result_s$optimal_lambda,
        Cd = C3,
        sigma2 = sigma2,
        max_iter = max_iter_nuclear,
        verbose = FALSE,
        warm_start = TRUE,
        use_newton = TRUE
      )

      M_debiased_list[[s]] <- M10_refined_list[[s]] + result_nuclear_s$Delta
    }

    # Average debiased estimates
    M_debiased_avg <- Reduce("+", M_debiased_list) / num_sources
    M_trans <- best_rank_r(M_debiased_avg, r = r)

    return(list(
      M_trans = M_trans,
      M_debiased_avg = M_debiased_avg,
      M_debiased_list = M_debiased_list,
      M10_refined_list = M10_refined_list,
      M_refined_list = M_refined_list,
      method = "DA",
      num_sources = num_sources
    ))

  } else {
    stop("Method must be either 'AD' or 'DA'")
  }
}

# ============================================================================
# Part 11: Source Detection Functions
# ============================================================================

#' Identify potential sources based on rank comparison using IC criterion
#'
#' @param X_sources List of source data matrices (may contain missing values)
#' @param X0 Target data matrix (may contain missing values)
#' @param r_max Maximum number of factors to consider (default: 10)
#' @param ic_type IC criterion type: "IC1" or "IC2" (default: "IC1")
#' @param data_type Type of data: "continuous", "count", or "binary"
#' @param C CJMLE projection constant (if NULL, auto-calculated)
#' @param max_iter Maximum CJMLE iterations (default: 30)
#' @param verbose Print progress information (default: TRUE)
#' @return List with positive_potential_sources, negative_sources, r_target, r_sources
#'
#' @examples
#' \donttest{
#' # Generate Poisson data
#' set.seed(2025)
#'
#' # Generate 5 sources with different ranks
#' n1 <- 100; p1 <- 100
#' source_list <- list()
#'
#' # Sources 1-2: rank 2 (same as target)
#' r_s <- 2
#' F_s <- matrix(runif(n1 * r_s, min = -2, max = 2), n1, r_s)
#' B_s <- matrix(runif(p1 * r_s, min = -2, max = 2), p1, r_s)
#' M_s <- F_s %*% t(B_s)
#' for (s in 1:2) {
#'   X_s <- matrix(rpois(n1 * p1, exp(M_s)), n1, p1)
#'
#'   # Add 10% missing values
#'   n_missing <- floor(n1 * p1 * 0.1)
#'   missing_idx <- sample(n1 * p1, n_missing)
#'   X_s[missing_idx] <- NA
#'
#'   source_list[[s]] <- X_s
#' }
#'
#' # Sources 3-5: rank 3 (different from target)
#' for (s in 3:5) {
#'   r_s_nega <- 3
#'   F_s_nega <- matrix(runif(n1 * r_s_nega, min = -2, max = 2), n1, r_s_nega)
#'   B_s_nega <- matrix(runif(p1 * r_s_nega, min = -2, max = 2), p1, r_s_nega)
#'   M_s_nega <- F_s_nega %*% t(B_s_nega)
#'   X_s_nega <- matrix(rpois(n1 * p1, exp(M_s_nega)), n1, p1)
#'
#'   n_missing <- floor(n1 * p1 * 0.1)
#'   missing_idx <- sample(n1 * p1, n_missing)
#'   X_s_nega[missing_idx] <- NA
#'
#'   source_list[[s]] <- X_s_nega
#' }
#'
#' # Target data: rank 2
#' n0 <- 50; p0 <- 50; r_target <- 2
#' M_target <- M_s[1:n0, 1:p0]
#' X_target <- matrix(rpois(n0 * p0, exp(M_target)), n0, p0)
#'
#' # Identify potential sources
#' result <- source_potential(
#'   X_sources = source_list,
#'   X0 = X_target,
#'   r_max = 5,
#'   ic_type = "IC1",
#'   data_type = "count",
#'   verbose = TRUE
#' )
#'
#' print(result$positive_potential_sources)  # Should be c(1, 2)
#' print(result$negative_sources)            # Should be c(3, 4, 5)
#' print(result$r_target)                    # Should be 2
#' print(result$r_sources)                   # Should be c(2, 2, 3, 3, 3)
#' }
#' @export
source_potential <- function(X_sources, X0, r_max = 10,
                             ic_type = "IC1", data_type = "count",
                             C = NULL, max_iter = 30, verbose = TRUE) {

  sourcenumber <- length(X_sources)

  if (verbose) {
    cat("\n===== Source Potential Detection (Rank-Based) =====\n")
    cat(sprintf("Source number %d\n", sourcenumber))
    cat(sprintf("IC type %s\n", ic_type))
  }

  # Step 1: Estimate target rank
  if (verbose) cat("\nStep 1 Target data \n")

  ic_result_target <- ic_criterion(X0, r_max = r_max, ic_type = ic_type,
                                   data_type = data_type, C = C,
                                   max_iter = max_iter, verbose = verbose)
  r_target <- ic_result_target$r_hat

  if (verbose) {
    cat(sprintf("Target estimated r = %d\n", r_target))
  }

  # Step 2: Estimate each source's rank
  if (verbose) cat("\nStep 2 Source data \n")

  r_sources <- numeric(sourcenumber)
  ic_results_sources <- list()

  for (i in 1:sourcenumber) {
    if (verbose) cat(sprintf("\n  Source %d:\n", i))

    ic_result_source <- ic_criterion(X_sources[[i]], r_max = r_max,
                                     ic_type = ic_type, data_type = data_type,
                                     C = C, max_iter = max_iter,
                                     verbose = verbose)
    r_sources[i] <- ic_result_source$r_hat
    ic_results_sources[[i]] <- ic_result_source

    if (verbose) {
      cat(sprintf("  Source %d estimated r = %d\n", i, r_sources[i]))
    }
  }

  # Step 3: Classify based on rank
  if (verbose) cat("\nStep 3 Classify \n")

  positive_potential_sources <- which(r_sources == r_target)
  negative_sources <- which(r_sources != r_target)

  if (verbose) {
    cat(sprintf("\nTarget number r = %d\n", r_target))
    cat("Source number \n")
    for (i in 1:sourcenumber) {
      cat(sprintf("  Source %d: r = %d %s\n",
                  i, r_sources[i],
                  ifelse(r_sources[i] == r_target, "(positive potential)", "(negative)")))
    }
    cat(sprintf("\nPositive Potential sources (r = %d): [%s]\n",
                r_target, paste(positive_potential_sources, collapse = ", ")))
    if (length(negative_sources) > 0) {
      cat(sprintf("Negative sources (r not equal %d): [%s]\n",
                  r_target, paste(negative_sources, collapse = ", ")))
    } else {
      cat("Negative sources: []\n")
    }
  }

  return(list(
    positive_potential_sources = positive_potential_sources,
    negative_sources = negative_sources,
    r_target = r_target,
    r_sources = r_sources,
    ic_type = ic_type,
    ic_result_target = ic_result_target,
    ic_results_sources = ic_results_sources
  ))
}

#' Detect positive and negative transfer sources using ratio method
#'
#' @param X_sources List of source data matrices (may contain missing values)
#' @param X0 Target data matrix (complete)
#' @param r Number of factors
#' @param C CJMLE projection constant
#' @param C2 Refinement projection constant
#' @param data_type Type of data: "continuous", "count", or "binary"
#' @param c_penalty Penalty coefficient (default: 0.1)
#' @param verbose Print progress information (default: TRUE)
#' @return List with positive_sources, negative_sources, and diagnostic info
#' @export
source_detection <- function(X_sources, X0, r, C, C2, data_type = "count",
                             c_penalty = 0.1, verbose = TRUE) {
  sourcenumber <- length(X_sources)
  n0 <- nrow(X0)
  p0 <- ncol(X0)

  if (verbose) {
    cat("\n===== Source Detection (Ratio Method) =====\n")
    cat(sprintf("Source number  %d\n", sourcenumber))
    cat(sprintf("Target dim %d x %d\n", n0, p0))
    cat(sprintf("Data type %s\n", data_type))
  }

  # Step 1: Estimate each source (CJMLE + Refine)

  M_source_estimates <- list()
  source_missing_rates <- numeric(sourcenumber)

  for (i in 1:sourcenumber) {
    X_current <- X_sources[[i]]
    n1 <- nrow(X_current)
    p1 <- ncol(X_current)

    # Create Omega matrix
    Omega_current <- matrix(1, nrow = n1, ncol = p1)
    Omega_current[is.na(X_current)] <- 0
    source_missing_rates[i] <- 1 - mean(Omega_current)

    if (verbose) {
      cat(sprintf("  Source %d: missing rate = %.2f%%\n", i, source_missing_rates[i] * 100))
    }

    X_current[is.na(X_current)] <- 0

    # Initialize
    if (data_type == "count") {
      init_data <- log(X_current + 1)
    } else {
      init_data <- X_current
    }

    svd_init <- svd(init_data)
    d_sqrt <- sqrt(abs(svd_init$d[1:r]))
    F_init <- svd_init$u[, 1:r, drop = FALSE] %*% diag(d_sqrt, nrow = r)
    B_init <- svd_init$v[, 1:r, drop = FALSE] %*% diag(d_sqrt, nrow = r)

    # CJMLE
    result_cjmle <- CJMLE(F_init, B_init, r, X_current,
                          Omega_current, C, data_type, method = "newton")
    M_hat <- result_cjmle$M

    # Refine
    result_refine <- refi.nosp(M_hat, r, X_current,
                               Omega_current, C2, data_type)
    M_refined <- result_refine$M

    # Extract target region
    M_source_estimates[[i]] <- M_refined[1:n0, 1:p0]
  }

  # Step 2: Compute Naive estimate on target

  if (data_type == "count") {
    init_data0 <- log(X0 + 1)
  } else {
    init_data0 <- X0
  }

  svd_init0 <- svd(init_data0)
  d0_sqrt <- sqrt(abs(svd_init0$d[1:r]))
  F0_init <- svd_init0$u[, 1:r, drop = FALSE] %*% diag(d0_sqrt, nrow = r)
  B0_init <- svd_init0$v[, 1:r, drop = FALSE] %*% diag(d0_sqrt, nrow = r)

  result_naive <- CJMLE(F0_init, B0_init, r, X0,
                        matrix(1, nrow = n0, ncol = p0), C, data_type,
                        method = "newton")
  M_naive <- result_naive$M

  # Step 3: Compute log-likelihoods

  L0 <- lik(M_naive, X0, matrix(1, nrow = n0, ncol = p0), data_type)

  L_sources <- numeric(sourcenumber)
  for (i in 1:sourcenumber) {
    L_sources[i] <- lik(M_source_estimates[[i]], X0,
                        matrix(1, nrow = n0, ncol = p0), data_type)
  }

  if (verbose) {
    cat(sprintf("  L0 (Naive): %.4f\n", L0))
    cat("  Li (Sources):\n")
    for (i in 1:sourcenumber) {
      cat(sprintf("    L%d: %.4f\n", i, L_sources[i]))
    }
  }

  # Step 4: Compute criterion

  penalty <- c_penalty * max(n0, p0)
  criterion_values <- abs(L_sources - L0) + penalty

  if (verbose) {
    cat(sprintf(" Penalty: %.4f\n", penalty))
    cat(" Value:\n")
    for (i in 1:sourcenumber) {
      cat(sprintf("    Source %d: abs(%.4f - %.4f) + %.4f = %.4f\n",
                  i, L_sources[i], L0, penalty, criterion_values[i]))
    }
  }

  # Step 5: Sort and find gap using ratio


  sorted_indices <- order(criterion_values)
  sorted_values <- criterion_values[sorted_indices]

  if (verbose) {
    cat("Sorted source:\n")
    for (i in 1:sourcenumber) {
      cat(sprintf(" Located %d: Source %d, value =%.4f, missing rate=%.2f%%\n",
                  i, sorted_indices[i], sorted_values[i],
                  source_missing_rates[sorted_indices[i]] * 100))
    }
  }

  # Compute ratios between adjacent values
  if (sourcenumber > 1) {
    ratios <- numeric(sourcenumber - 1)
    for (i in 1:(sourcenumber - 1)) {
      if (sorted_values[i] < 1e-10) {
        ratios[i] <- sorted_values[i + 1] / 1e-10
      } else {
        ratios[i] <- sorted_values[i + 1] / sorted_values[i]
      }
    }

    # Find maximum gap
    gap_position <- which.max(ratios)

    if (verbose) {
      cat(" Ratio:\n")
      for (i in 1:(sourcenumber - 1)) {
        cat(sprintf(" Ratio %d (Source %d / Source %d): %.4f %s\n",
                    i, sorted_indices[i+1], sorted_indices[i], ratios[i],
                    if(i == gap_position) " max gap" else ""))
      }
      cat(sprintf(" Max gap location: %d\n", gap_position))
    }

    # Classify positive and negative
    positive_sources <- sorted_indices[1:gap_position]
    negative_sources <- sorted_indices[(gap_position + 1):sourcenumber]

  } else {
    gap_position <- NA
    positive_sources <- 1
    negative_sources <- integer(0)
  }

  if (verbose) {
    cat(sprintf("\n=====Results=====\n"))
    cat(sprintf("Positive transfer sources: [%s]\n",
                paste(positive_sources, collapse = ", ")))
    if (length(negative_sources) > 0) {
      cat(sprintf("Negative transfer sources: [%s]\n",
                  paste(negative_sources, collapse = ", ")))
    } else {
      cat("Negative transfer sources: []\n")
    }
  }

  return(list(
    positive_sources = positive_sources,
    negative_sources = negative_sources,
    loglik_diffs = L_sources - L0,
    criterion_values = criterion_values,
    sorted_indices = sorted_indices,
    sorted_values = sorted_values,
    gap_position = gap_position,
    L0 = L0,
    L_sources = L_sources,
    source_missing_rates = source_missing_rates
  ))
}

# ============================================================================
# Part 12: Information Criterion for Rank Selection
# ============================================================================

#' Information criterion (IC1/IC2) for selecting number of factors
#'
#' @param X Data matrix (may contain missing values coded as NA)
#' @param r_max Maximum number of factors to consider (default: 10)
#' @param ic_type IC criterion type: "IC1" or "IC2" (default: "IC1")
#' @param data_type Type of data: "continuous", "count", or "binary"
#' @param C CJMLE projection constant (if NULL, auto-calculated)
#' @param max_iter Maximum CJMLE iterations (default: 30)
#' @param verbose Print progress information (default: FALSE)
#' @return List with r_hat (optimal rank), ic_values, loglik_values
#'
#' @examples
#' # Generate Poisson data with known rank
#' set.seed(2025)
#' n <- 100; p <- 100; r_true <- 2
#'
#' # Generate true factors
#' F_true <- matrix(runif(n * r_true, min = -2, max = 2), n, r_true)
#' B_true <- matrix(runif(p * r_true, min = -2, max = 2), p, r_true)
#' M_true <- F_true %*% t(B_true)
#'
#' # Generate Poisson observations
#' lambda <- exp(M_true)
#' X <- matrix(rpois(n * p, as.vector(lambda)), n, p)
#'
#' # Add 10% missing values
#' n_missing <- floor(n * p * 0.1)
#' missing_idx <- sample(n * p, n_missing)
#' X[missing_idx] <- NA
#'
#' # Use IC1 to select rank
#' result_IC1 <- ic_criterion(
#'   X = X,
#'   r_max = 6,
#'   ic_type = "IC1",
#'   data_type = "count",
#'   verbose = TRUE
#' )
#'
#' print(paste("True rank:", r_true))
#' print(paste("Estimated rank (IC1):", result_IC1$r_hat))
#'
#' # Use IC2 to select rank
#' result_IC2 <- ic_criterion(
#'   X = X,
#'   r_max = 6,
#'   ic_type = "IC2",
#'   data_type = "count",
#'   verbose = TRUE
#' )
#'
#' @export
ic_criterion <- function(X, r_max = 10, ic_type = c("IC1", "IC2"),
                         data_type = "count", C = NULL,
                         max_iter = 30, verbose = FALSE) {

  ic_type <- match.arg(ic_type)

  n <- nrow(X)
  p <- ncol(X)

  # Create Omega matrix for missing values
  Omega <- matrix(1, nrow = n, ncol = p)
  Omega[is.na(X)] <- 0
  X[is.na(X)] <- 0

  if (verbose) {
    cat(sprintf("\n===== IC =====\n"))
    cat(sprintf("Date dim %d x %d\n", n, p))
    cat(sprintf("Missing rate %.2f%%\n", (1 - mean(Omega)) * 100))
    cat(sprintf("IC type %s\n", ic_type))
    cat(sprintf("Candidate r from 1 to %d\n", r_max))
  }

  # Ensure r_max does not exceed min(n,p)
  r_max <- min(r_max, n, p)

  # If C not provided, auto-calculate
  if (is.null(C)) {
    C <- ceiling(sqrt(r_max))
  }

  # Store results
  ic_values <- numeric(r_max)
  loglik_values <- numeric(r_max)

  # Compute IC for each candidate r
  for (r in 1:r_max) {
    if (verbose) cat(sprintf("Calculate r=%d...\n", r))

    # Initialize
    if (data_type == "count") {
      init_data <- log(X + 1)
    } else {
      init_data <- X
    }

    svd_init <- svd(init_data)
    d_sqrt <- sqrt(abs(svd_init$d[1:r]))
    F_init <- svd_init$u[, 1:r, drop = FALSE] %*% diag(d_sqrt, nrow = r)
    B_init <- svd_init$v[, 1:r, drop = FALSE] %*% diag(d_sqrt, nrow = r)

    # CJMLE estimation
    tryCatch({
      result_cjmle <- CJMLE(F_init, B_init, r, X, Omega, C,
                            data_type, method = "newton", max_iter = max_iter)
      M_hat <- result_cjmle$M

      # Compute log-likelihood
      loglik <- lik(M_hat, X, Omega, data_type)
      loglik_values[r] <- loglik

      # Compute penalty
      if (ic_type == "IC1") {
        # IC1: max(n,p) * log(min(n,p))
        penalty <- max(n, p) * log(min(n, p))
      } else {
        # IC2: (n+p) * log(n*p/(n+p))
        penalty <- (n + p) * log(n * p / (n + p))
      }

      # IC = -2*loglik + r*penalty
      ic_values[r] <- -2 * loglik + r * penalty

      if (verbose) {
        cat(sprintf("    r=%d: loglik=%.4f, penalty=%.4f, IC=%.4f\n",
                    r, loglik, r * penalty, ic_values[r]))
      }
    }, error = function(e) {
      if (verbose) cat(sprintf("    r=%d: fail \n", r))
      ic_values[r] <- Inf
      loglik_values[r] <- -Inf
    })
  }

  # Select r with minimum IC
  r_hat <- which.min(ic_values)

  if (verbose) {
    cat(sprintf("\n Estimated number: r = %d\n", r_hat))
    cat(sprintf("IC value: %.4f\n", ic_values[r_hat]))
    cat(sprintf("Log-like value: %.4f\n", loglik_values[r_hat]))
  }

  return(list(
    r_hat = r_hat,
    ic_values = ic_values,
    loglik_values = loglik_values,
    ic_type = ic_type
  ))
}

# ============================================================================
# End of transGFM Function Library
# ============================================================================
