#' Convert cd/m2 to dB lambert
#'
#' @param x_b A number representing a luminance in cd/m2.
#' @returns A number representing the respective luminance in dB Lambert
#' relative to 10^-10 Lambert.
#' @export
db_lambert <- function(x_b) {
  10 * log10(x_b * pi / (10^(-6)))
}

#' Convert dB lambert to cd/m2
#'
#' @param db_b A number representing a luminance in dB Lambert.
#' @returns A number representing a luminance in cd/m2.
#' @export
db_inv_lambert <- function(db_b) {
  10^(db_b / 10) * 10^(-6) / pi
}

#' Convert sound pressure in Pascal to dB SPL re 20 muPa
#'
#' @param x_l A number representing a sound pressure in Pascal.
#' @returns A number representing the respective sound in dB SPL relative
#' to 20 muPa.
#' @export
db_spl <- function(x_l) {
  20 * log10(x_l / (2 * 10^(-5)))
}

#' Convert dB SPL to Pascal
#'
#' @param db_l A number representing a sound in dB SPL.
#' @returns A number representing the sound pressure in Pascal.
#' @export
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
#' @export
psi <- function(x, alpha, beta) {
  return(alpha * x^beta)
  }

#' Inverse of psychophysical function phi(x) = alpha * x^beta
#'
#' @param x A number representing a perceived intentsity (e.g.
#' brightness or loudness)
#' @param alpha A number.
#' @param beta A number.
#' @returns A number representing a physical stimulus intentsity (e.g.
#' luminance in cd/m2 or sound pressure in Pascal).
#' @export
psi_inv <- function(x, alpha, beta) {
  return((1 / alpha * x)^(1 / beta))
}

#' Cognitive weighting function W(p)
#' @export
weigh_fun <- function(p, w_1, w = .6) (w_1 * p^w)

#' Calculate sum of internal references
#' @export
const_fun <- function(w_p, rho_std, db_inv_std, alpha_std, beta_std,
                           rho_tgt, db_inv_tgt, alpha_tgt, beta_tgt) {
  w_p * psi(db_inv_std(rho_std), alpha_std, beta_std) -
        psi(db_inv_tgt(rho_tgt), alpha_tgt, beta_tgt)
}

#' Calculate global psychophysics model prediction of the physical intensity
#' produced in a magnitude production task with one production ratio
#'
#' @param standard_intensity A number or vector representing the physical
#' intensity/ies of the standard stimulus / stimuli (in dB).
#' @param alpha_std A number for the parameter of the psychophysical function
#' of the standard.
#' @param alpha_tgt A number for the parameter of the psychophysical function
#' of the target.
#' @param beta_std A number for the parameter of the psychophysical function
#' of the standard.
#' @param beta_tgt A number for the parameter of the psychophysical function
#' of the target.
#' @param w_p A number representing the weighted production ratios W(p).
#' @param rho_std A number for the internal reference parameter pertaining to
#' the standard dimension.
#' @param rho_tgt A number for the internal reference parameter pertaining to
#' the target dimension.
#' @param const A number representing the sum of rho's in the restricted model
#'        for matching / p=1 only.
#' @param task Either "bright_loud" or "loud_bright" indicating the task.
#' @returns A number or vector representing the physical intensity predicted by
#' the global psychophysics model.
#' @export
# TODO this should also work for a matrix of params?
gpm <- function(standard_intensity,
                alpha_std, alpha_tgt,
                beta_std, beta_tgt,
                w_p,
                rho_std = NULL, rho_tgt = NULL, const = NULL,
                task = c("bright_loud", "loud_bright")) {

  # Setzt alpha_b reliabel auf 1 aka das alpha, das nicht übergeben wird
  if (is.null(alpha_std) || length(alpha_std) == 0) alpha_std <- 1
  if (is.null(alpha_tgt) || length(alpha_tgt) == 0) alpha_tgt <- 1
  
  if (task == "bright_loud") {
    db_inv_std <- db_inv_lambert
    db_inv_tgt <- db_inv_spl
    db_tgt <- db_spl
  } else {
    db_inv_std <- db_inv_spl
    db_inv_tgt <- db_inv_lambert
    db_tgt <- db_lambert
  }
  if(!(is.null(rho_std) & is.null(rho_tgt))) {
    const <- const_fun(w_p, rho_std, db_inv_std, alpha_std, beta_std,
                           rho_tgt, db_inv_tgt, alpha_tgt, beta_tgt)
  }else if(is.null(rho_std)){
    stop(glue::glue("The internal reference of the standard stimuli in {task} is NULL. Please provide a valid input for this parameter"))
  }else{
    stop(glue::glue("The internal reference of the target stimuli in {task} is NULL. Please provide a valid input for this parameter"))
  }
  (w_p * psi(db_inv_std(standard_intensity), alpha_std, beta_std) - const) |>
    psi_inv(alpha_tgt, beta_tgt) |>
    db_tgt() |>
    stats::setNames(rep("x_p", length(standard_intensity)))
}

#' Calculate global psychophysics model prediction 
#'
#' @description 
#' This function extends the global psychophysics model function (gpm) by 
#' allowing calculation of model predictions for a specific production ratio 
#' (p). It handles the weighting calculation internally and provides appropriate 
#' naming of results based on the production ratio used.
#' 
#' @param standard_intensity A number or vector representing the physical
#'        intensity/ies of the standard stimulus/stimuli (in dB).
#' @param alpha_std A number for the parameter of the psychophysical function
#'        of the standard.
#' @param alpha_tgt A number for the parameter of the psychophysical function
#'        of the target.
#' @param beta_std A number for the parameter of the psychophysical function
#'        of the standard.
#' @param beta_tgt A number for the parameter of the psychophysical function
#'        of the target.
#' @param p An integer (1, 2, or 3) representing the production ratio to be used.
#' @param w_1 A number representing the weight parameter for the weighting 
#'            function.
#' @param w A number representing the exponent in the weighting function. 
#'          Default is 0.6.
#' @param w_p A direct specification of the weighted production ratio. 
#'        If w_1 and w are given, this will be ignored.
#' @param rho_std A number for the internal reference parameter pertaining to
#'        the standard dimension.
#' @param rho_tgt A number for the internal reference parameter pertaining to
#'        the target dimension.
#' @param const A number representing the sum of rho's in the restricted model
#'        for matching / p=1 only.
#' @param task Either "bright_loud" or "loud_bright" indicating the task.
#'
#' @returns A number or vector representing the physical intensity predicted by
#'          the global psychophysics model, with names reflecting the production 
#'          ratio used.
#'
#' @examples
#' # Calculate prediction for brightness-to-loudness with production ratio p=2
#' gpm_multiple_p(standard_intensity = 70,
#'                alpha_std = 1, alpha_tgt = 0.5,
#'                beta_std = 0.33, beta_tgt = 0.67,
#'                rho_std = 20, rho_tgt = 20,
#'                p = 2, w_1 = 1, w = 0.6,
#'                task = "bright_loud")
#' @export
gpm_multiple_p <- function(standard_intensity,
                        alpha_std, alpha_tgt,
                        beta_std, beta_tgt,
                        p, w_1 = NULL, w = NULL, w_p = NULL,
                        rho_std = NULL, rho_tgt = NULL, const = NULL,
                        task = c("bright_loud", "loud_bright")) {
  
  # Validate production factor p
  if (!(p %in% 1:3))
    stop("Production factor p not valid, must be 1, 2, or 3")
  
  # Calculate w_p based on p, w_1, and w
  if (!is.null(w_1) && !is.null(w)) {
    w_p <- weigh_fun(p, w_1, w)
  } else if (!is.null(w_p)) {
    w_p <- w_p
  } else {
    stop("Either w_p or both w_1 and w must be provided")
  }
  
  # Call the original gpm function with calculated w_p
  result <- gpm(standard_intensity = standard_intensity,
                alpha_std = alpha_std, alpha_tgt = alpha_tgt,
                beta_std = beta_std, beta_tgt = beta_tgt,
                w_p = w_p, 
                rho_std = rho_std, rho_tgt = rho_tgt, 
                const = const,
                task = task)
  
  names(result) <- rep(paste0("x_", p), length(result))
  
  return(result)
}
