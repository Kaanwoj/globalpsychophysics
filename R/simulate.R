#' @importFrom stats rbeta rlnorm rnorm runif
#' @examples
#' param <- create_param_set()
#' param <- create_param_set(restriction = "role-independent")
#' param <- create_param_set(restriction = "const")
#' @export
create_param_set <- function(n = 1, 
                             method = "existing", 
                             restriction = "no",
                             calc_omega_p = FALSE) {
  # Match the argument with valid choices and get the first matching value
  method <- match.arg(method, choices = c("existing", "generate"))
  restriction <- match.arg(restriction, choices = c("no", "role-independent", "const"))

  if (method == "existing") {
    out <- data.frame(alpha_b = 1, alpha_l = 15.48,
                      beta_b = 0.4, beta_l = 0.26)
    if (restriction == "no") {
      out$rho_btol <- 38.3
      out$rho_lfromb <- 20
      out$rho_ltob <- 32.4
      out$rho_bfroml <- 71.8
    }
    if (restriction == "role-independent") {
      out$rho_b <- 78.2
      out$rho_l <- 51.4
    }
    if (restriction == "const") {
      out$const_bl <- -1.622
      out$const_lb <- 0.0803
    }
    if (calc_omega_p) {
      out$omega_p <- 0.8
    } else {
      out$omega <- 0.6
      out$omega_1 <- 0.8
    }
  } else {
    out <- data.frame(alpha_b = 1,
                      alpha_l = rlnorm(n, 12, 3),
                      beta_b = rbeta(n, 3, 6),
                      beta_l = rbeta(n, 3, 6))
    if (restriction == "no") {
      out$rho_btol <- runif(n, 50, 90)
      out$rho_lfromb <- runif(n, 0, 80)
      out$rho_ltob <- runif(n, 0, 80)
      out$rho_bfroml <- runif(n, 50, 90)
    }
    if (restriction == "role-independent") {
      out$rho_b <- runif(n, 50, 90)
      out$rho_l <- runif(n, 0, 80)
    }
    if (restriction == "const") {
      out$const_bl <- runif(n, 100, 180)
      out$const_lb <- runif(, 0, 160)
    }
    if (calc_omega_p) {
      out$omega_p <- rnorm(n, 1, .3)
    } else {
      out$omega <- rnorm(n, 0.6, .05)
      out$omega_1 <- rnorm(n, 0.8, .1)
    }
  }
  out
}

#' Predict magnitude productions from global psychophysics model
#'
#' @description
#' A wrapper function that extracts parameters from a parameter set and calls
#' the appropriate global psychophysics model function to predict magnitude productions.
#'
#' @param standards A vector of numbers representing the physical intensities 
#'                  of standard stimuli (in dB).
#' @param task A character string, either "bright_loud" or "loud_bright" 
#'             indicating the task.
#' @param param A named vector of numbers representing the model parameters.
#' @param p An integer (default: 1) representing the production ratio to be used.
#' @returns A vector of predicted physical intensities.
#' @examples
#' param <- create_param_set()
#' y_bright <- predict_gpm(c(34, 43, 52, 70), "loud_bright", param, p = 2)
#' y_loud <- predict_gpm(c(69, 73, 77, 85), "bright_loud", param)
#' @export
predict_gpm <- function(standards, task, param, p = 1) {
  if ("omega_1" %in% names(param)) {
    omega_1 <- param$omega_1
  } else {
    omega_1 <- NULL 
  }
  if ("omega_p" %in% names(param)) {
    omega_p <- param$omega_p
  } else {
    omega_p <- NULL
  }
  if ("omega" %in% names(param)) {
    omega <- param$omega
  } else {
    omega <- NULL 
  }
  
  if((is.null(param$rho_bfroml) | is.null(param$rho_lfromb) | is.null(param$rho_ltob) | is.null(param$rho_btol))){
    if(task == "bright_loud"){
      rho_std <- param$rho_b
      rho_tgt <- param$rho_l
      alpha_std <- param$alpha_b
      alpha_tgt <- param$alpha_l
      beta_std <- param$beta_b
      beta_tgt <- param$beta_l
    }else{
      rho_std <- param$rho_l
      rho_tgt <- param$rho_b
      alpha_std <- param$alpha_l
      alpha_tgt <- param$alpha_b
      beta_std <- param$beta_l
      beta_tgt <- param$beta_b
    }
  }else{
    if(task == "bright_loud"){
      rho_std <- param$rho_btol
      rho_tgt <- param$rho_lfromb
      alpha_std <- param$alpha_b
      alpha_tgt <- param$alpha_l
      beta_std <- param$beta_b
      beta_tgt <- param$beta_l
    }else{
      rho_std <- param$rho_ltob
      rho_tgt <- param$rho_bfroml
      alpha_std <- param$alpha_l
      alpha_tgt <- param$alpha_b
      beta_std <- param$beta_l
      beta_tgt <- param$beta_b
    }
  }

  gpm_multiple_p(standard_intensity = standards,
              alpha_std = alpha_std, alpha_tgt = alpha_tgt,
              beta_std = beta_std, beta_tgt = beta_tgt,
              p = as.numeric(p), 
              omega_1 = omega_1, omega_p = omega_p, omega = omega,
              task = task,
              rho_std = rho_std, rho_tgt = rho_tgt)
}

#' Simulate a magnitude production data set for both tasks and one or more
#' production ratios
#' 
#' @description
#' This function generates simulated magnitude production data by predicting responses 
#' using the global psychophysics model and adding random noise.
#'
#' @param ntrials A number of trials per condition.
#' @param cond A data frame with columns:
#'        - task: factor with values "loud_bright" or "bright_loud"
#'        - std: numeric values representing standard intensities
#'        - sigma: standard deviation for adding noise to predictions
#'        - p: production ratio values
#' @param param A named vector of numbers representing the model parameters:
#'                alpha_b, alpha_l, beta_b, beta_l,
#'                omega_p or omega_1 and omega (if omega_p should be calculated) 
#'                rho_btol, rho_lfromb, rho_ltob, rho_bfroml if role-dependece 
#'                or rho_b and rho_l if role-independece is assumed
#' @returns A data frame with ntrials rows per condition, containing columns:
#'          task, std, sigma, p, mu (predicted value), and tgt (noisy response).
#' @examples
#' Single production factor example 
#' cond <- data.frame(
#'                    task = rep(c("loud_bright", "bright_loud"), each = 6),
#'                    std = c(25, 34, 43, 52, 61, 70,  # loud std
#'                            69, 73, 77, 81, 85, 89), # bright std
#'                    sigma = rep(c(bright = 3,  loud = 6), each = 6),
#'                    p = c(1)
#' )
#' param <- create_param_set()
#' simulate_gpm(20, cond, param)
#' 
#' Multiple production factors example
#' param <- create_param_set()
#' 
#' loud_standards <- c(25, 34, 43, 52, 61, 70)
#' bright_standards <- c(69, 73, 77, 81, 85, 89)
#' 
#' cond_loud <- expand.grid(task = "loud_bright",
#'                          std = loud_standards,
#'                          p = c(1, 2, 3))
#'                          
#' cond_bright <- expand.grid(task = "bright_loud",
#'                            std = bright_standards,
#'                            p = c(1, 2, 3))
#'                            
#' cond_multi_p <- rbind(cond_loud, cond_bright)
#' cond_multi_p$sigma <- ifelse(cond_multi_p$task == "bright_loud", 3, 3)
#' 
#' simulated_data <- simulate_gpm(200, cond_multi_p, param)
#' @export
simulate_gpm <- function(ntrials, cond, param) {
  cond$task <- factor(cond$task, levels = unique(cond$task))
  
  # Group by unique combination of task and p
  cond_split <- split(cond, list(cond$task, cond$p))
  
  predictions <- lapply(cond_split, function(x) {
    # Each group has the same task and p value
    task_val <- as.character(x$task[1])
    p_val <- x$p[1]
    
    # Generate predictions for this group
    preds <- predict_gpm(
      standards = x$std, 
      task = task_val,
      param = param,
      p = p_val
    )
    
    # Return a data frame with row indices and predictions
    data.frame(
      row_idx = match(x$std, cond$std[cond$task == task_val & cond$p == p_val]),
      mu = preds
    )
  })
  
  # Combine predictions and assign to original conditions by index
  mu_values <- numeric(nrow(cond))
  for (i in seq_along(predictions)) {
    name_parts <- strsplit(names(cond_split)[i], "\\.")[[1]]
    task_name <- name_parts[1]
    p_value <- as.numeric(name_parts[2])
    
    group_rows <- which(cond$task == task_name & cond$p == p_value)
    pred_df <- predictions[[i]]
    mu_values[group_rows] <- pred_df$mu
  }
  cond$mu <- mu_values
  
  # Create output data frame with replicated rows for trials
  out <- cond[rep(seq_len(nrow(cond)), each = ntrials), ]
  
  # Check for invalid mu values
  invalid_mu <- is.na(out$mu) | is.infinite(out$mu)
  if(any(invalid_mu)){
    warning(paste("Invalid mu values detected for", sum(invalid_mu), 
                  "rows. Check your predict_gpm function outputs."))
    out$mu[invalid_mu] <- 50
  }
  
  # Check for invalid sigma values
  invalid_sigma <- is.na(out$sigma) | out$sigma < 0
  if(any(invalid_sigma)){
    warning(paste("Invalid sigma values detected for", sum(invalid_sigma), 
                  "rows. Using default value."))
    out$sigma[invalid_sigma] <- 1  
  }
  
  out$tgt <- rnorm(nrow(out), mean = out$mu, sd = out$sigma)
  row.names(out) <- NULL
  out
}

