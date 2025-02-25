test_that("gpm_bright_loud_1", {
  expect_equal(
               gpm(standard_intensity = 70,
                   alpha_std = 1, alpha_tgt = 15.48,
                   beta_std = 0.4, beta_tgt = 0.26,
                   w_p = 0.8,
                   task = "bright_loud",
                   rho_std = NULL, rho_tgt = NULL, const = -1.51) |>
                 round(2),
               c(x_p = 36.63))
})

test_that("gpm_vector", {
  expect_equal(
               gpm(standard_intensity = c(70, 70),
                   alpha_std = 1, alpha_tgt = 15.48,
                   beta_std = 0.4, beta_tgt = 0.26,
                   w_p = 0.8,
                   task = "bright_loud",
                   rho_std = NULL, rho_tgt = NULL, const = -1.51) |>
                 round(2),
               c(x_p = 36.63, x_p = 36.63))
})
