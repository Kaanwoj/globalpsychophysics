library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#' Create a list of data information to be given to rstan::stan()
#'
#' @param data A data frame containing columns standard and match in dB, and
#' task direction (bright_loud or loud_bright).
#' @param ntrials number of trials per standard.
#' @returns
#' @export
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
  
  # sort by p and then by std -------------------------------------------------
  if ("p" %in% names(data)) {   # multiple production factors
    data <- sort_by(data, ~ p + std)
  } else {
    data <- sort_by(data, ~ std)
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

  # check task ----------------------------------------------------------------
  if (!"task" %in% names(data)) {
    if ("type" %in% names(data)) {
      data$task <- data$type
    } else {
      stop("column 'task' or 'type' is missing in data")
    }
  }

  # split by task -------------------------------------------------------------
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
    datlist$nstd_p_lb <- length(datlist$std_bl) * datlist$np
    datlist$nstd_p_bl <- length(datlist$std_lb) * datlist$np
   #datlist$tgt_lb <- dat_lb$tgt
   #datlist$tgt_bl <- dat_bl$tgt
    datlist$tgt <- data$tgt
    datlist$sig_lb <- aggregate(tgt ~ std + p, dat_lb, sd)$tgt
    datlist$sig_bl <- aggregate(tgt ~ std + p, dat_bl, sd)$tgt
   #datlist$sig <- aggregate(tgt ~ std + p + task, data, sd)$tgt
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

#' Build a Stan data list for the global psychophysics model
#'
#' @description
#' Constructs a named list suitable for passing directly to \code{rstan::stan()}
#' for the three-modality global psychophysics model (\code{gpm_generic.stan}).
#' Handles modality encoding, condition indexing, sigma indexing, and prior
#' hyperparameters. Accepts either real experimental data or an
#'  conditions data frame.
#'
#' @param dat A data frame of experimental data (already converted to dB).
#'   Must contain columns for standard intensity, target modality, standard
#'   modality, production factor, and produced intensity. See Details for
#'   accepted column names. If \code{NULL}, \code{cond} must be provided.
#' @param cond A data frame of experimental conditions for prior predictive
#'   checks, e.g. built with \code{expand.grid}. No produced intensity column
#'   required - \code{tgt} will be filled with zeros and \code{onlyprior}
#'   forced to \code{1}. If \code{NULL}, \code{dat} must be provided.
#' @param prior_params A named list of prior hyperparameters to pass to Stan,
#'   e.g. \code{list(alpha_l_logmu = 4, alpha_l_logsigma = 0.3)}. Any
#'   parameters not supplied are filled from \code{default_priorset.csv}. 
#'   Defaults to an empty list, in which case
#'   all priors come from the defaults.
#' @param onlyprior Integer, \code{0} or \code{1}. If \code{1}, the likelihood
#'   is skipped in Stan (prior predictive mode) and \code{tgt} is set to zeros
#'   regardless of whether \code{dat} is provided. Automatically set to
#'   \code{1} when only \code{cond} is given. Default \code{0}.
#' @param sigidx An optional integer vector of length \code{nrow(dat)} or
#'   \code{nrow(cond)} mapping each trial to a sigma parameter index. If
#'   \code{NULL} (default), the function checks for a \code{sigidx} column
#'   in \code{dat}/\code{cond}; if absent, all trials are assigned index
#'   \code{1}.
#' @param tgt_col Character. Name of the column in \code{dat} containing the
#'   produced intensity values in dB. Default \code{"matchdb"}.
#'
#' @details
#' \strong{Modality encoding:} Modalities are encoded as integers:
#' \code{1} = loudness (auditory), \code{2} = brightness (visual),
#' \code{3} = vibration strength (tactile). The function accepts either a
#' \code{task} column with strings like \code{"loud_bright"} or separate
#' \code{standard_modality} and \code{target_modality} columns with strings
#' \code{"auditory"}/\code{"loud"}, \code{"visual"}/\code{"bright"},
#' \code{"tactile"}/\code{"strong"}.
#'
#' \strong{Accepted column names:}
#' \itemize{
#'   \item Standard intensity: \code{std}, \code{standard}, or \code{std_db}
#'   \item Production factor: \code{p} or \code{prod_factor}
#'   \item Sigma index: \code{sigidx} (optional, see \code{sigidx} parameter)
#' }
#'
#' \strong{Condition indexing:} Multiple trials sharing the same combination
#' of \code{std_modality}, \code{tgt_modality}, \code{p}, and \code{std} are
#' assigned the same condition index \code{idx}, which maps trials to the
#' corresponding \code{mu} in Stan.
#'
#' @return A named list ready to pass to \code{rstan::stan()} as the
#'   \code{data} argument.
#'
#' @examples
#' # Prior predictive check with expand.grid conditions
#' cond <- rbind(
#'   expand.grid(task = "loud_bright", std = c(25, 34, 43, 52, 61, 70), p = 1:3),
#'   expand.grid(task = "bright_loud", std = c(61, 63, 67, 72, 76, 81), p = 1:3)
#' )
#' datlist <- make_datlist(cond = cond)
#'
#' # Real data with custom prior set
#' prior_sets <- read.table("priorsets.csv", header = TRUE, sep = ";")
#' param <- setNames(as.list(prior_sets$set_4), prior_sets$parname)
#' datlist <- make_datlist(dat = df, prior_params = param, tgt_col = "match_db")
#'
#' # Real data with sigma indexing via column in dat
#' dat$sigidx <- case_when(
#'   dat$task == "loud_bright" & dat$std <= 25 ~ 1,
#'   dat$task == "loud_bright" & dat$std >= 61 ~ 3,
#'   dat$task == "loud_bright"                 ~ 2,
#'   dat$task == "bright_loud" & dat$std <= 63 ~ 4,
#'   dat$task == "bright_loud" & dat$std >= 81 ~ 6,
#'   dat$task == "bright_loud"                 ~ 5
#' )
#' datlist <- make_datlist(dat = dat, prior_params = param, tgt_col = "match_db")
#'
#' @export

make_datlist_generic <- function(
    dat = NULL,          # real data frame from data files
    cond = NULL,         # expand.grid-style conditions frame
    prior_params = list(),        # named list/vector of prior parameters
    onlyprior = 0,
    sigidx = NULL,
    tgt_col = "matchdb"
) {
  
  library(dplyr)
  
  # --- 1. Modality string -> integer mapping --------------------------------
  modality_to_int <- function(x) {
    case_when(
      x == "auditory" | x == "loud"    ~ 1,
      x == "visual"   | x == "bright"  ~ 2,
      x == "tactile"  | x == "strong"  ~ 3
    )
  }
  
  # --- 2. Resolve which data frame we're working with ----------------------
  # cond path: prior predictive, no real tgt values
  # dat path:  real data, has tgt_col
  if (!is.null(cond) && is.null(dat)) {
    df <- cond
    tgt_vals <- numeric(nrow(df))  # dummy zeros
    onlyprior <- 1                 # force prior predictive if only cond given
    message("No data given, prior predicitve checks will be performed")
  } else if (!is.null(dat) && is.null(cond)) {
    df <- dat
    if(onlyprior == 1){
      tgt_vals <- numeric(nrow(df))  # dummy zeros
    }else{
      tgt_vals <- dat[[tgt_col]]
    }
  } else {
    stop("Provide either dat or cond, not both or neither.")
  }
  
  # --- 3. Resolve std_modality and tgt_modality ----------------------------
  # Accept either task string ("loud_bright") or separate modality columns
  if ("task" %in% colnames(df)) {
    # parse from task e.g. "loud_bright" -> std = "loud", tgt = "bright"
    parts <- strsplit(as.character(df$task), "_")
    std_modality <- modality_to_int(sapply(parts, '[', 1))
    tgt_modality <- modality_to_int(sapply(parts, '[', 2))
  } else if (all(c("standard_modality", "target_modality") %in% colnames(df))) {
    std_modality <- modality_to_int(df$standard_modality)
    tgt_modality <- modality_to_int(df$target_modality)
  } else {
    stop("Provide either a 'task' column or both 'standard_modality' and 
         'target_modality' columns.")
  }
  
  # --- 4. Resolve std column -----------------------------------------------
  # accept either "std" or "standard" as column name
  if ("std" %in% colnames(df)) {
    std_vals <- df$std
  } else if ("standard" %in% colnames(df)) {
    std_vals <- df$standard
  } else if ("std_db" %in% colnames(df)) {
    std_vals <- df$std_db
  } else {
    stop("No standard intensity column found. Expected 'std' or 'standard'.")
  }
  
  # --- 5. Resolve p column -------------------------------------------------
  if ("p" %in% colnames(df)) {
    p_vals <- df$p
  } else if ("prod_factor" %in% colnames(df)) {
    p_vals <- df$prod_factor
  } else {
    stop("No production factor column found. Expected 'p' or 'prod_factor'.")
  }
  
  # --- 6. Build condition index --------------------------------------------
  # multiple trials per unique (std_modality, tgt_modality, p, std) share idx
  cond_key <- data.frame(std_modality, tgt_modality, p = p_vals, std = std_vals)
  idx <- cond_key %>%
    group_by(std_modality, tgt_modality, p, std) %>%
    mutate(idx = cur_group_id()) %>% # assigns unique index
    pull(idx)
  
  # unique conditions for ncond
  cond_unique <- cond_key %>%
    group_by(std_modality, tgt_modality, p, std) %>%
    summarise(.groups = "drop") %>%
    arrange(std_modality, tgt_modality, p, std)
  
  # --- 7. Build sigidx -----------------------------------------------------
  if (!is.null(sigidx)) {
    sigidx_vals <- sigidx                     
  } else if ("sigidx" %in% colnames(df)) {
    sigidx_vals <- df$sigidx                   
  } else {
    # Default: same sigma index for everything 
    sigidx_vals <- rep(1, nrow(df))          
  }
  
  # --- Fill in unused priors with default ----------------------------------
  default_prior_params <- read.table("default_priors.csv", 
                                     header = TRUE, sep = ";") %>%
    (\(x) setNames(as.list(x$value), x$parname))()
  
  resolved_params <- modifyList(default_prior_params, prior_params)
  
  # --- 8. Assemble datlist -------------------------------------------------
  datlist <- c(
    list(
      ntotal       = nrow(df),
      ncond        = nrow(cond_unique),
      nsig         = max(sigidx_vals),
      std_modality = cond_unique$std_modality,
      tgt_modality = cond_unique$tgt_modality,
      p            = cond_unique$p,
      std          = cond_unique$std,
      tgt          = tgt_vals,
      idx          = idx,
      sigidx       = sigidx_vals,
      onlyprior    = onlyprior
    ),
    resolved_params  # prior parameters appended directly
  )
  
  return(datlist)
}


#' Estimate parameters of global psychophysics model
#'
#' @param data A data frame containing the standard and target intensity in dB,
#' and task direction (bright_loud or loud_bright). Must be sorted loud_bright
#' first and then bright_loud.
#' @param ntrials number of trials per standard.
#' @returns
#' @examples
#' data(matching)
#' datlist <- make_datlist(data = matching[matching$id == 2, ], 48)
#' m <- estimate(data = matching[matching$id == 2, ], 48, references = "constant")
#' print(m, pars = c("alpha_l", "beta_l", "beta_b", "omega1", "omega", "const_lb",
#'                   "const_bl"), probs = c(.025, .975))
#' 
#' @export
estimate <- function(data, ntrials,
                     references = c("dependent", "independent", "constant")) {
  references <- match.arg(references)
  # TODO: args for rstan::stan()

  datlist <- make_datlist(data, ntrials)

  if ("p" %in% names(datlist)) {    # multiple p model
    print("fitting model with multiple production factors and role-dependent
      internal references")
      # FIXME: hard coded path
    model <- "R/stan/gpm.stan"
  } else {                          # p=1 models
    if (references == "constant") {
      print("fitting model with constant sum of internal references")
      model <- "R/stan/gpm_p1_const.stan"
    } else {
      print("fitting model with role-dependent internal references")
      model <- "R/stan/gpm_p1.stan"
    }
  }

  stan(file = model, data = datlist, chains = 4, warmup = 1000, iter = 10000,
       control = list(max_treedepth = 15))
}
