if(!require(globalpsychophysics)){
  devtools::install_github("https://github.com/Kaanwoj/globalpsychophysics.git")
}

param <- parametersets |> subset(set == "set_3") |>
  _[, c("alpha_l", "beta_b", "beta_l", "rho_btol", "rho_bfroml", "rho_ltob",
        "rho_lfromb", "omega_1", "omega")]

param$rho_btol <- 70
ntrials <- 3
p <- c(1, 2, 3)
loud_standards <- c(25)
bright_standards <- c(69)
cond_basic <- rbind(expand.grid(task = "loud_bright", std = loud_standards, p = p ),
              expand.grid(task = "bright_loud", std = bright_standards, p = p))
cond_basic$sigma <- 3

dat_basic <- simulate_gpm(ntrials, cond_basic, param, replace_invalid = FALSE)

test_that("simulate_check_dimensions", {
  expect_equal(
    dim(dat_basic), c(18, 6))
})

test_that("simulate_sound_vals", {
  expect_equal(
    gpm(
      standard_intensity=bright_standards,
      alpha_std=NULL,
      alpha_tgt=param$alpha_l,
      beta_std=param$beta_b,
      beta_tgt=param$beta_l,
      omega_p=weigh_fun(p=1, param$omega_1, param$omega),
      rho_std=param$rho_btol,
      rho_tgt=param$rho_lfromb,
      const=NULL,
      task = "bright_loud") |> unname(),
    dat_basic[dat_basic$task == "bright_loud" & dat_basic$p == 1, "mu"][1])
})

test_that("simulate_light_vals", {
  expect_equal(
    sapply(p, \(pi) {
      gpm(
        standard_intensity=loud_standards,
        alpha_std=param$alpha_l,
        alpha_tgt=NULL,
        beta_std=param$beta_l,
        beta_tgt=param$beta_b,
        omega_p=weigh_fun(p=pi, param$omega_1, param$omega),
        rho_std=param$rho_ltob,
        rho_tgt=param$rho_bfroml,
        const=NULL,
        task = "loud_bright")
      }) |> unname(),
    dat_basic[dat_basic$task == "loud_bright", "mu"][seq(1,
      ntrials*length(p)*length(loud_standards), by = ntrials)])
})
