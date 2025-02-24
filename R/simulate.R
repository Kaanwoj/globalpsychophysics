#' @importFrom stats rbeta rlnorm rnorm runif
create_param_set <- function(n = 1, method = "existing", restriction = "no") {
  if (method == "existing") {
    out <- data.frame(alpha_b = 1, alpha_l = 15.48,
             beta_b = 0.4, beta_l = 0.26,
             w_1 = 0.8,
             rho_btol = 60.5, rho_lfromb = 21.8,
             rho_ltob = 8.6, rho_bfroml = 60.5)
  } else {
    out <- data.frame(alpha_b = 1,
             alpha_l = rlnorm(n, 12, 3),
             beta_l = rbeta(n, 3, 6),
             beta_b = rbeta(n, 3, 6),
             w_1 = rnorm(n, 1, .3),
             rho_btol = runif(n, 50, 90), rho_lfromb = runif(n, 0, 80),
             rho_ltob = runif(n, 0, 80), rho_bfroml = runif(n, 50, 90))
  }
  out
}

#' Predict data for global psychophysics model
#'
#' @param task A character string, either "bright_loud" or "loud_bright"
#' indicating the task.
#' @param standards A vector of numbers.
#' @param param A named vector of numbers representing the model parameters.
#' @returns A vector of predicted values.
#' @examples
#' par <- create_param_set()
#' y_loud <- predict_gpm("bright_loud", c(60, 70, 77, 80), par)
#' y_bright <- predict_gpm("loud_bright", c(34, 45, 56, 69), par)
predict_gpm <- function(task, standards, param) {
  gpm("bright_loud", standards,
      alpha_std = param$alpha_b, alpha_tgt = param$alpha_l,
      beta_std = param$beta_b, beta_tgt = param$beta_l,
      w_p = param$w_1,
      rho_std = param$rho_btol, rho_tgt = param$rho_lfromb)
}
