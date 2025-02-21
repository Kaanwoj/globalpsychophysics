#' Cross-modal matching of brightness and loudness
#'
#' A data set of seven subjects perfoming matching tasks adjusting the
#' brightness of visual stimuli to match the perceived intensity of noise
#' sounds and vice versa. Six standard levels are used in each dimension and
#' per standard each subject made 48 matches. Each row is from one trial.
#'
#' @format ## `matching`
#' A data frame with 4,032 rows and 6 columns:
#' \describe{
#'   \item{id}{Subject id}
#'   \item{type}{Type of matching direction, i.e. whether sound (bright_loud)
#'   or light (loud_bright) is adjusted}
#'   \item{light}{Luminance in dB Lambert}
#'   \item{sound}{Sound pressure in dB SPL}
#'   \item{standard_db}{Intensity of standard stimulus in dB}
#'   \item{match_db}{Intensity of final adjustment of the target stimulus in dB}
#' }
"matching"
