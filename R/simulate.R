#' @importFrom stats rbeta rlnorm rnorm runif
#' @examples
#' param <- create_param_set()
#' param <- create_param_set(restriction = "role-independent")
#' param <- create_param_set(restriction = "const")
#' @export
create_param_set <- function(n = 1,
                             method = "existing",
                             restriction = "no",
                             calc_omega_p = FALSE,
                             task = c("bright_loud", "loud_bright",
                                      "strong_loud", "loud_strong",
                                      "bright_strong", "strong_bright")){
  
  method      <- match.arg(method, choices = c("existing", "generate"))
  restriction <- match.arg(restriction, choices = c("no", "role-independent", 
                                                    "const"))
  task        <- match.arg(task)
  
  standard_modality <- sub("_.*", "", task)
  target_modality   <- sub(".*_", "", task)
  modalities        <- unique(c(standard_modality, target_modality))
  
  if (method == "existing") {
    # TODO insert reasonable values for _t parameters
    out <- data.frame(alpha_b = 1, alpha_l = 15.48, alpha_t = 1,
                      beta_b  = 0.4, beta_l = 0.26, beta_t = 0.33)
    
    # Keep only alphas/betas relevant to the task
    alpha_cols <- paste0("alpha_", substring(modalities, 1, 1))
    beta_cols  <- paste0("beta_",  substring(modalities, 1, 1))
    out <- out[, c(alpha_cols, beta_cols), drop = FALSE]
    
    if (restriction == "no") {
      std <- substring(standard_modality, 1, 1)
      tgt <- substring(target_modality,   1, 1)
      out[[paste0("rho_", std, "to", tgt)]]   <- switch(task,
                                                        bright_loud   = 38.3,  
                                                        loud_bright  = 32.4,
                                                        strong_loud   = 35,    
                                                        loud_strong  = 32.4,
                                                        bright_strong = 38.3,  
                                                        strong_bright = 35)
      out[[paste0("rho_", tgt, "from", std)]] <- switch(task,
                                                        bright_loud   = 20,    
                                                        loud_bright  = 71.8,
                                                        strong_loud   = 50,    
                                                        loud_strong  = 71.8,
                                                        bright_strong = 20,    
                                                        strong_bright = 50)
    }
    if (restriction == "role-independent") {
      for (m in modalities) {
        out[[paste0("rho_", substring(m, 1, 1))]] <- switch(m,
                                                            bright = 78.2, 
                                                            loud = 51.4, 
                                                            strong = 60)
      }
    }
    if (restriction == "const") {
      std <- substring(standard_modality, 1, 1)
      tgt <- substring(target_modality,   1, 1)
      out[[paste0("const_", std, tgt)]] <- switch(task,
                                                  bright_loud   = -1.622, 
                                                  loud_bright  =  0.0803,
                                                  strong_loud   =  0.01,     
                                                  loud_strong  =  0.0803,
                                                  bright_strong =  -1.622,     
                                                  strong_bright = 0.01)
    }
    if (calc_omega_p) {
      out$omega_p <- 0.8
    } else {
      out$omega   <- 0.6
      out$omega_1 <- 0.8
    }
    
  } else {
    # Generate mode
    out <- data.frame(row.names = seq_len(n))
    
    # Add alpha/beta only for relevant modalities
    for (m in modalities) {
      initial <- substring(m, 1, 1)
      if("bright" %in% modalities){
        out[[paste0("alpha_", initial)]] <- switch(m,
                                                   bright = rep(1, n),
                                                   loud   = rlnorm(n, 12, 3),
                                                   strong = rlnorm(n, 12, 3))
      }else{
        out[[paste0("alpha_", initial)]] <- switch(m,
                                                   loud   = rlnorm(n, 12, 3),
                                                   strong = rep(1, n))
      }
      
      out[[paste0("beta_", initial)]] <- switch(m,
                                                bright = rbeta(n, 3, 6),
                                                loud   = rbeta(n, 3, 6),
                                                strong = rbeta(n, 3, 6))
    }
    
    std <- substring(standard_modality, 1, 1)
    tgt <- substring(target_modality,   1, 1)
    
    if (restriction == "no") {
      out[[paste0("rho_", std, "to", tgt)]]   <- switch(standard_modality,
                                                    bright = runif(n, 0, 80),
                                                    loud   = runif(n, 50, 90),
                                                    strong = runif(n, 50, 90))
      
      out[[paste0("rho_", tgt, "from", std)]] <- switch(target_modality,
                                                    bright = runif(n, 0, 80),
                                                    loud   = runif(n, 50, 90),
                                                    strong = runif(n, 50, 90))
      out[[paste0("rho_", tgt, "to", std)]]   <- switch(target_modality,
                                                    bright = runif(n, 0, 80),
                                                    loud   = runif(n, 50, 90),
                                                    strong = runif(n, 50, 90))
      out[[paste0("rho_", std, "from", tgt)]] <- switch(standard_modality,
                                                    bright = runif(n, 0, 80),
                                                    loud   = runif(n, 50, 90),
                                                    strong = runif(n, 50, 90))
    }
    if (restriction == "role-independent") {
      for (m in modalities) {
        out[[paste0("rho_", substring(m, 1, 1))]] <- switch(m,
                                                    bright = runif(n, 0, 80),
                                                    loud   = runif(n, 50, 90),
                                                    strong = runif(n, 50, 90))
      }
    }
    if (restriction == "const") {
      # If-loop only to preserve old runif range for lb
      if(paste0(std, tgt) == "lb"){
        out$const_lb <- runif(n, 0, 160)
      }else{
        out[[paste0("const_", std, tgt)]] <- runif(n, 100, 180)  
      }
    }
    if (calc_omega_p) {
      out$omega_p <- rnorm(n, 1,   .3)
    } else {
      out$omega   <- rnorm(n, 0.6, .05)
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
#' @param task String indicating the task, e.g "bright_loud", "loud_bright", 
#' "strong_loud", "bright_strong", etc.
#' @param param A named vector of numbers representing the model parameters.
#' @param p An integer (default: 1) representing the production ratio to be used.
#' @returns A vector of predicted physical intensities.
#' @examples
#' param <- create_param_set()
#' y_bright <- predict_gpm(c(34, 43, 52, 70), "loud_bright", param, p = 2)
#' y_loud <- predict_gpm(c(69, 73, 77, 85), "bright_loud", param)
#' @export
predict_gpm <- function(standards, task, param, p = 1) {
  
  std <- substring(sub("_.*", "", task), 1, 1)
  tgt <- substring(sub(".*_", "", task), 1, 1)
  
  # Extract omega parameters
  omega_1 <- param$omega_1 %||% NULL
  omega_p <- param$omega_p %||% NULL
  omega   <- param$omega   %||% NULL
  
  # Determine rho/const based on restriction type
  rho_std_col   <- paste0("rho_", std, "to", tgt)
  rho_tgt_col   <- paste0("rho_", tgt, "from", std)
  rho_std_ri    <- paste0("rho_", std)
  rho_tgt_ri    <- paste0("rho_", tgt)
  const_col     <- paste0("const_", std, tgt)
  
  if (!is.null(param[[rho_std_col]]) & !is.null(param[[rho_tgt_col]])) {
    # in case of role-dependence
    rho_std <- param[[rho_std_col]]
    rho_tgt <- param[[rho_tgt_col]]
    const   <- NULL
  } else if (!is.null(param[[rho_std_ri]]) & !is.null(param[[rho_tgt_ri]])) {
    # in case of role-independence
    rho_std <- param[[rho_std_ri]]
    rho_tgt <- param[[rho_tgt_ri]]
    const   <- NULL
  } else if (!is.null(param[[const_col]])) {
    rho_std <- NULL
    rho_tgt <- NULL
    const   <- param[[const_col]]
  } else {
    stop(glue::glue(
      "No valid parameters found for task '{task}'. Provide one of:
       - role-dependent: {rho_std_col} and {rho_tgt_col}
       - role-independent: {rho_std_ri} and {rho_tgt_ri}
       - constant sum: {const_col}"
    ))
  }
  
  # Extract alpha/beta by modality
  alpha_std <- param[[paste0("alpha_", std)]]
  alpha_tgt <- param[[paste0("alpha_", tgt)]]
  beta_std  <- param[[paste0("beta_",  std)]]
  beta_tgt  <- param[[paste0("beta_",  tgt)]]
  
  gpm_multiple_p(standard_intensity = standards,
                 alpha_std = alpha_std, alpha_tgt = alpha_tgt,
                 beta_std  = beta_std,  beta_tgt  = beta_tgt,
                 p = as.numeric(p),
                 omega_1 = omega_1, omega = omega, omega_p = omega_p,
                 rho_std = rho_std, rho_tgt = rho_tgt,
                 const   = const,
                 task    = task)
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
#'        - task: factor indicating the task, e.g "bright_loud", "loud_bright", 
#'                "strong_loud", "bright_strong", etc.
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
  
  out <- cond[rep(seq_len(nrow(cond)), each = ntrials), ]
  
  # Define sensible fallback mu per target modality
  fallback_mu <- list(
    loud   = 15,   # dB SPL
    bright = 52,   # dB Lambert
    strong = 5     # dB displacement 
  )
  
  # Check for invalid mu values
  invalid_mu <- is.na(out$mu) | is.infinite(out$mu)
  if (any(invalid_mu)) {
    message("--- BEFORE ---")
    print(head(out[invalid_mu, ]))
    
    # Get target modality for each row and look up fallback
    tgt_modality <- sub(".*_", "", as.character(out$task))
    fb <- unlist(fallback_mu[tgt_modality])
    
    warning(paste(sum(invalid_mu), "invalid mu values detected and replaced",
                  "with modality-specific fallback values.",
                  "Check your predict_gpm function outputs."))
    
    out$mu <- ifelse(invalid_mu, fb, out$mu)
    
    message("--- AFTER ---")
    print(head(out[invalid_mu, ]))
  }
  
  # Check for invalid sigma values
  invalid_sigma <- is.na(out$sigma) | out$sigma < 0
  if (any(invalid_sigma)) {
    warning(paste(sum(invalid_sigma), "invalid sigma values detected and set to 1."))
    out$sigma[invalid_sigma] <- 1
  }
  
  out$tgt <- rnorm(nrow(out), mean = out$mu, sd = out$sigma)
  row.names(out) <- NULL
  out
}

