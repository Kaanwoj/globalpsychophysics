
# Setup -----
make_conditions <- function() {
  loud_standards   <- c(25, 34, 43, 52, 61, 70)
  strong_standards <- c(5, 10, 15, 20, 25, 30)
  
  cond_loud <- expand.grid(
    task = "loud_strong",
    std  = loud_standards,
    p    = c(1, 2, 3)
  )
  cond_bright <- expand.grid(
    task = "strong_loud",
    std  = strong_standards,
    p    = c(1, 2, 3)
  )
  
  cond <- rbind(cond_loud, cond_bright)
  cond$sigma <- 3
  cond
}

make_params <- function() {
  create_param_set(
    n           = 1,
    method      = "generate",
    restriction = "no",
    task        = "strong_loud"
  )
}

# create_param_set()  -----
test_that("create_param_set() returns a non-empty object", {
  param <- make_params()
  expect_true(!is.null(param))
  expect_true(length(param) > 0)
})

test_that("create_param_set() contains at least one alpha, two betas, one omega, 
          two rhos", {
  param <- make_params()
  nms   <- names(param)
  
  count_prefix <- function(prefix) sum(grepl(paste0("^", prefix, "_"), nms))
  
  expect_gte(count_prefix("alpha"), 1)
  expect_gte(count_prefix("beta"),  2)
  expect_gte(count_prefix("omega"), 1)
  expect_gte(count_prefix("rho"),   2)
})


# Condition setup  -----
test_that("condition data frame has correct structure", {
  cond <- make_conditions()
  expect_s3_class(cond, "data.frame")
  expect_named(cond, c("task", "std", "p", "sigma"), ignore.order = TRUE)
})

# simulate_gpm()  -----
test_that("simulate_gpm() returns a data frame", {
  param <- make_params()
  cond  <- make_conditions()
  sim   <- simulate_gpm(200, cond, param)
  expect_s3_class(sim, "data.frame")
})

test_that("simulate_gpm() returns expected columns", {
  param <- make_params()
  cond  <- make_conditions()
  sim   <- simulate_gpm(200, cond, param)
  expect_true(all(c("task", "std", "p", "mu", "tgt") %in% names(sim)))
})

test_that("simulate_gpm() is reproducible with set.seed()", {
  param <- make_params()
  cond  <- make_conditions()
  set.seed(42); sim1 <- simulate_gpm(10, cond, param)
  set.seed(42); sim2 <- simulate_gpm(10, cond, param)
  expect_equal(sim1$tgt, sim2$tgt)
})


