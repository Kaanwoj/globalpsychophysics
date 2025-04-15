
plot_matching <- function(param, limbright = c(10, 100), limloud = c(10, 100)) {
  splmin <- limloud[1]
  splmax <- limloud[2]
  lambertmin <- limbright[1]
  lambertmax <- limbright[2]
  plot(1, type = "n", xlim = c(60, 100), ylim = c(20, 100),
       xlab = "luminance [dB Lambert]", ylab = "sound pressure level [dB SPL]")
  points(lambertmin:lambertmax,
         gpm(lambertmin:lambertmax, param$alpha_b, param$alpha_l, param$beta_b,
             param$beta_l, param$w_p, rho_std = param$rho_btol,
             rho_tgt = param$rho_lfromb,
             task = "bright_to_loud"),
         type = "l", col = "steelblue4", lwd = 2)
  lines(splmin:splmax ~ gpm(splmin:splmax, param$alpha_l, param$alpha_b,
                            param$beta_l, param$beta_b, param$w_p,
                            rho_std = param$rho_ltob,
                            rho_tgt = param$rho_bfroml,
                            task = "loud_to_bright"),
        type = "l", col = "red3", lwd = 2)
  abline(h = c(25, 80), v = c(65, 89), lty = 2, col = "gray")
  abline(h = c(param$rho_slb, param$rho_cbl), v = c(param$rho_sbl, param$rho_clb), lty = 3, col = "gray")
  legend("bottomright", c(expression(brightness %->% loudness),
                  expression(loudness %->% brightness)),
         col = c("steelblue4", "red3"), pch = c(15, 19), lty = 1)
}

#' Plot target intensity values by standard intensity across different production factors
#' 
#' @description
#' Creates a faceted plot comparing the relationship between standard intensity and 
#' target intensity across different production factors (p) and directions. 
#' Each facet represents a different task, with lines of different colors 
#' representing different production factor values.
#' 
#' @param simulated_data A data frame containing simulation results with at least the 
#'        following columns:
#'        \itemize{
#'          \item \code{task}: Factor indicating different tasks
#'          \item \code{std}: Standard intensity values
#'          \item \code{p}: Production factor values 
#'          \item \code{tgt}: Target intensity values
#'        }
#' 
#' @import ggplot2
#' 
#' @examples
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
#' 
#' # Generate the plot
#' plot_multiple_p(simulated_data)
#' 
#' @export
plot_multiple_p <- function(simulated_data) {
  # Aggregate data to get mean and sd
  summary_data <- aggregate(
    tgt ~ task + std + p, 
    data = simulated_data,
    FUN = function(x) c(mean = mean(x), sd = sd(x))
  )
  
  # Reshape for plotting
  summary_data$mean <- summary_data$tgt[,1]
  summary_data$sd <- summary_data$tgt[,2]
  
  custom_colors <- c("#951A36", "#164467", "#76B3AD")
  # Plot comparing simulated data
  ggplot(summary_data, aes(x = std, y = mean, color = factor(p))) +
    geom_line() +
    geom_point() +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2) +
    facet_wrap(~ task, scales = "free") +
    scale_color_manual(values = custom_colors) +
    labs(x = "Standard Intensity",
         y = "Target Intensity",
         color = "Production\nFactor") +
    theme_minimal()
}
