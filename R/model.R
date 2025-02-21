#' Convert cd/m2 to dB lambert
#'
#' @param x_b A number representing a luminance in cd/m2.
#' @returns A number representing the respective luminance in dB Lambert
#' relative to 10^-10 Lambert.
db_lambert <- function(x_b) {
  10 * log10(x_b * pi / (10^(-6)))
}

#' Convert dB lambert to cd/m2
#'
#' @param db_b A number representing a luminance in dB Lambert.
#' @returns A number representing a luminance in cd/m2.
db_inv_lambert <- function(db_b) {
  10^(db_b / 10) * 10^(-6) / pi
}

#' Convert sound pressure in Pascal to dB SPL re 20 muPa
#'
#' @param x_l A number representing a sound pressure in Pascal.
#' @returns A number representing the respective sound in dB SPL relative
#' to 20 muPa.
db_spl <- function(x_l) {
  20 * log10(x_l / (2 * 10^(-5)))
}

#' Convert dB SPL to Pascal
#'
#' @param db_l A number representing a sound in dB SPL.
#' @returns A number representing the sound pressure in Pascal.
db_inv_spl <- function(db_l) {
  10^(db_l / 20) * 2 * 10^-5
}

#' Psychophysical function phi(x) = alpha * x^beta
#'
#' @param x A number representing a physical stimulus intentsity (e.g.
#' luminance in cd/m2 or sound pressure in Pascal).
#' @param alpha A number.
#' @param beta A number.
#' @returns A number representing a perceived intentsity (e.g.
#' brightness or loudness)
psi <- function(x, alpha, beta) alpha * x^beta

#' Inverse of psychophysical function phi(x) = alpha * x^beta
#'
#' @param x A number representing a perceived intentsity (e.g.
#' brightness or loudness)
#' @param alpha A number.
#' @param beta A number.
#' @returns A number representing a physical stimulus intentsity (e.g.
#' luminance in cd/m2 or sound pressure in Pascal).
psi_inv <- function(x, alpha, beta) (1 / alpha * x)^(1 / beta)

#' Calculate global psychophysics model prediction of the phsyical intensity
#' produced in a magnitude production task
#'
#' @param task Either "bright_loud" or "loud_bright" indicating the task.
#' @param standard_intensity A number or vector representing the physical
#' intensity or intensities of the standard stimulus / stimuli (in dB).
#' @param A numeric vector with the parameter values for ...
#' @param db_inv_std A function corresponding to the standard dimension.
#' @param db_inv_tgt A function corresponding to the standard dimension.
#' @param db_tgt A function corresponding to the standard dimension.
#' @returns A number or vector representing the physical intensity predicted by
#' the global psychophysics model.
# TODO this should also work for a matrix of params?
gpm <- function(task, standard_intensity,
                alpha_std, alpha_tgt,
                beta_std, beta_tgt,
                w_p,
                rho_std = NULL, rho_tgt = NULL, const = NULL) {
  # set alpha_b to 1, and corresponding db functions
  if (task == "bright_loud") {
    alpha_std <- 1
    db_inv_std <- db_inv_lambert
    db_inv_tgt <- db_inv_spl
    db_tgt <- db_spl
  } else {
    alpha_tgt <- 1
    db_inv_std <- db_inv_spl
    db_inv_tgt <- db_inv_lambert
    db_tgt <- db_lambert
  }
  if (is.null(rho_std)) { # restricted model for matching / p=1 only
    # TODO reduce code duplication here
    (w_p * psi(db_inv_std(standard_intensity), alpha_std, beta_std) - const) |>
      psi_inv(alpha_tgt, beta_tgt) |>
      db_tgt() |>
      stats::setNames(rep("x_p", length(standard_intensity)))
  } else {
    (w_p * psi(db_inv_std(standard_intensity), alpha_std, beta_std) -
       psi(db_inv_std(rho_std), alpha_std, beta_std) +
       psi(db_inv_tgt(rho_std), alpha_tgt, beta_tgt)) |>
      psi_inv(alpha_tgt, beta_tgt) |>
      db_tgt() |>
      stats::setNames(rep("x_p", length(standard_intensity)))
  }
}
