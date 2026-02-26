


path <- "D:/Google Drive/promotion/projekte/globalpsychophysics"
devtools::load_all(path = path)  # in Rpackage globalpsychophysics folder
# https://github.com/Kaanwoj/globalpsychophysics

library(ggplot2)

# param <- list(
#   alpha_l = 15,     
#   beta_b = 0.4,      
#   beta_l = .2,   
#   omega_p = 0.7,
#   omega_1 = 0.7,
#   omega = 2,
#   rho_btol = 0.3,    
#   rho_lfromb = 0.3  
# )

param <- create_param_set(n = 1, method = "generate", restriction = "no",
                          task = "strong_loud")

# Multiple production factors example
loud_standards <- c(25, 34, 43, 52, 61, 70)
strong_standards <- c(5, 10, 15, 20, 25, 30)

cond_loud <- expand.grid(
  task = "loud_strong",
  std = loud_standards,
  p = c(1, 2, 3)
)

cond_bright <- expand.grid(
  task = "strong_loud",
  std = strong_standards,
  p = c(1, 2, 3)
)


cond_multi_p <- rbind(cond_loud, cond_bright)
cond_multi_p$sigma <- ifelse(cond_multi_p$task == "loud_strong", 3, 3)

simulated_data <- simulate_gpm(200, cond_multi_p, param)

# plot raw data ---------------------------------------------------------------
ggplot(simulated_data, aes(x = std, y = tgt, color = p, shape = task)) +
  geom_point() +
  labs(x = "Standard",
       y = "Target") +
  theme_minimal() +
  theme(legend.position="none")


# Aggregated plot -------------------------------------------------------------
summary_data <- aggregate(
  tgt ~ task + std + p, 
  data = simulated_data,
  FUN = function(x) c(mean = mean(x), sd = sd(x))
)
# Reshape for plotting
summary_data$mean <- summary_data$tgt[,1]
summary_data$sd <- summary_data$tgt[,2]

# Plot comparing p values 
ggplot(summary_data, aes(x = std, y = mean, color = factor(p))) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2) +
  facet_wrap(~ task, scales = "free") +
  labs(x = "Standard Intensity",
       y = "Target Intensity",
       color = "Production\nFactor") +
  theme_minimal()


# Summary of mus --------------------------------------------------------------
dat_summary <- unique(simulated_data[, c("task", "std", "p", "mu")])
head(dat_summary)

