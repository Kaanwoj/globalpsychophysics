#' @importFrom stats rbeta rlnorm rnorm runif
#' @examples
#' param <- create_param_set()
#' param <- create_param_set(restriction = "role-independent")
#' param <- create_param_set(restriction = "const")
create_param_set <- function(n = 1, method = c("existing", "generate"),
                             restriction = c("no", "role-independent", "const")) {
  method <- method[1]
  restriction <- restriction[1] # FIXME: is this how it should be dealt with?
  if (method == "existing") {
    out <- data.frame(alpha_b = 1, alpha_l = 15.48,
                      beta_b = 0.4, beta_l = 0.26,
                      w_p = 0.8)
    if (restriction == "no") {
      out$rho_btol <- 38.3
      out$rho_lfromb <- 20
      out$rho_ltob <- 32.4
      out$rho_bfroml <- 71.8
    }
    if (restriction == "role-independent") {
      out$rho_bfroml <- out$rho_btol <- 78.2
      out$rho_lfromb <- out$rho_ltob <- 51.4
    }
    if (restriction == "const") {
      out$const_bl <- -1.622
      out$const_lb <- 0.0803
    }
  } else {
    out <- data.frame(alpha_b = 1,
                      alpha_l = rlnorm(n, 12, 3),
                      beta_l = rbeta(n, 3, 6),
                      beta_b = rbeta(n, 3, 6),
                      w_p = rnorm(n, 1, .3),
                      rho_btol = runif(n, 50, 90), rho_lfromb = runif(n, 0, 80),
                      rho_ltob = runif(n, 0, 80), rho_bfroml = runif(n, 50, 90))
  }
  out
}

#' Predict magnitude productions from global psychophysics model
#'
#' @param task A character string, either "bright_loud" or "loud_bright"
#' indicating the task.
#' @param standards A vector of numbers.
#' @param param A named vector of numbers representing the model parameters.
#' @returns A vector of predicted values.
#' @examples
#' param <- create_param_set()
#' y_bright <- predict_gpm(c(34, 43, 52, 70), "loud_bright", param)
#' y_loud <- predict_gpm(c(69, 73, 77, 85), "bright_loud", param)
predict_gpm <- function(standards, task, param, p = 1) {
  gpm(standards,
      alpha_std = param$alpha_b, alpha_tgt = param$alpha_l,
      beta_std = param$beta_b, beta_tgt = param$beta_l,
      w_p = param$w_p,
      task,
      rho_std = param$rho_btol, rho_tgt = param$rho_lfromb)
}

#' Simulate a magnitude production data set for both tasks and one or more
#' production ratios
#'
#' @param ntrials A number of trials per condition.
#' @param cond A data frame with columns task, std, sigma, and p.
#'        sigma: number representing the standard deviation for loudness and
#'        brightness productions.
#' @param param A named vector of numbers representing the model parameters.
#' @returns A data frame with the simulated magnitude productions for ntrials
#' per condition.
#' @examples
#' cond <- data.frame(
#'                    task = rep(c("loud_bright", "bright_loud"), each = 6),
#'                    std = c(25, 34, 43, 52, 61, 70,  # loud std
#'                            69, 73, 77, 81, 85, 89), # bright std
#'                    sigma = rep(c(bright = 3,  loud = 6), each = 6),
#'                    p = c(1)
#' )
#' param <- create_param_set()
#' simulate_gpm(20, cond, param)
# TODO make it work for other p's than 1
simulate_gpm <- function(ntrials, cond, param) {
  cond$task <- factor(cond$task, levels = unique(cond$task))
  cond$mu <- split(cond, cond$task) |>
    lapply(\(x) {predict_gpm(x$std, as.character(x$task)[1], param, p)}) |>
    unlist() |> unname()
  out <- cond[rep(seq_len(nrow(cond)), each = ntrials), ]
  out$tgt <- rnorm(nrow(out), mean = out$mu, sd = out$sigma)
  row.names(out) <- NULL
  out
}
