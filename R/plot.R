
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
