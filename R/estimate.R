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

  if (!"task" %in% names(data)) {
    if ("type" %in% names(data)) {
      data$task <- data$type
    } else {
      stop("column 'task' or 'type' is missing in data")
    }
  }

  dat_lb <- data |> subset(task == "loud_bright")
  std_lb <- unique(dat_lb$std)
  dat_bl <- data |> subset(task == "bright_loud")
  std_bl <- unique(dat_bl$std)

  # TODO: option for sigma to be a free parameter

  datlist <- list(
    ntotal = nrow(data),
    ntotal_lb = nrow(dat_lb),
    ntotal_bl = nrow(dat_bl),
    nstd_lb = std_lb |> length(),
    nstd_bl = std_bl |> length(),
    std_lb = std_lb,
    std_bl = std_bl
  )

  if ("p" %in% names(data)) {   # multiple production factors
    datlist$np <- length(unique(data$p))
    datlist$p <- unique(data$p)
    datlist$std_p_lb <- as.numeric(interaction(dat_lb$p, dat_lb$std,
                                     lex.order = TRUE))
    datlist$std_p_bl <- as.numeric(interaction(dat_bl$p, dat_bl$std,
                                     lex.order = TRUE))
    datlist$tgt <- data$tgt
    datlist$sig <- aggregate(tgt ~ std + p + task, data, sd)$tgt
  } else {                      # matching
    datlist$std_lb_idx <- as.numeric(as.factor(dat_lb$std))
    datlist$std_bl_idx <- as.numeric(as.factor(dat_bl$std))
    datlist$tgt_lb <- dat_lb$tgt
    datlist$tgt_bl <- dat_bl$tgt
    datlist$sig_lb <- aggregate(tgt ~ std, dat_lb, sd)$tgt
    datlist$sig_bl <- aggregate(tgt ~ std, dat_bl, sd)$tgt
  }

  datlist
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
  }
  stan(file = model, data = datlist, chains = 4, warmup = 1000, iter = 2000)
}
