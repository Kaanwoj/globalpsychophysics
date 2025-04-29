library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#' Create a list of data information to be given to rstan::stan()
#'
#' @param data A data frame containing columns standard and match in dB, and
#' task direction (bright_loud or loud_bright).
#' @param ntrials number of trials per standard.
#' @returns
make_datlist <- function(data, ntrials) {

  # check standard ------------------------------------------------------------
  if (!"std" %in% names(data)) {
    if ("standard_db" %in% names(data)) {
      data$std <- data$standard_db
    } else if ("standard" %in% names(data)) {
      data$std <- data$standard
    } else {
      stop("column 'standard', 'standard_db' or 'std' is missing in data")
    }
  }

  # check match ---------------------------------------------------------------
  if (!"tgt" %in% names(data)) {
    if ("match_db" %in% names(data)) {
      data$tgt <- data$match_db
    } else if ("match" %in% names(data)) {
      data$tgt <- data$match
    } else {
      stop("column 'match', 'match_db' or 'tgt' is missing in data")
    }
  }

  dat_lb <- data |> subset(type == "loud_bright")
  std_lb <- unique(dat_lb$std)
  dat_bl <- data |> subset(type == "bright_loud")
  std_bl <- unique(dat_bl$std)

  # TODO: option for sigma to be a free parameter

  list(ntotal_lb = nrow(dat_lb),
       ntotal_bl = nrow(dat_bl),
       ntrials = ntrials,
       nx_lb = std_lb |> length(),
       nx_bl = std_bl |> length(),
       x_lb = std_lb,
       x_bl = std_bl,
       x_lb_idx = as.numeric(as.factor(dat_lb$std)),
       x_bl_idx = as.numeric(as.factor(dat_bl$std)),
       y_lb = dat_lb$tgt,
       y_bl = dat_bl$tgt,
       sig_lb = aggregate(tgt ~ std, dat_lb, sd)$tgt,
       sig_bl = aggregate(tgt ~ std, dat_bl, sd)$tgt
  )
}

#' Estimate parameters of global psychophysics model
#'
#' @param data A data frame containing columns standard and match in dB, and
#' task direction (bright_loud or loud_bright).
#' @param ntrials number of trials per standard.
#' @returns
#' @examples
#' data(matching)
#' datlist <- make_datlist(data = matching[matching$id == 2, ], 48)
#' m <- estimate(data = matching[matching$id == 2, ], 48, references = "constant")
#' print(m, pars = c("alpha_l", "beta_l", "beta_b", "w_1", "const_lb",
#'                   "const_bl"), probs = c(.025, .975))
#' 
estimate <- function(data, ntrials,
                     references = c("dependent", "independent", "constant")) {
  # TODO: args for rstan::stan()
  
  datlist <- make_datlist(data, ntrials)

  if (references == "constant") {
    print("fitting model with constant sum of internal references")
    model <- "R/stan/gpm_p1_const.stan"
    # FIXME: no hard coded path
  } else {
    print("fitting model with role-dependent internal references")
    model <- "R/stan/gpm_p1.stan"
    datlist$alpha_b <- 1
  }
  stan(file = model, data = datlist, chains = 4, warmup = 1000, iter = 2000)
}
