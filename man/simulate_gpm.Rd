% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate.R
\name{simulate_gpm}
\alias{simulate_gpm}
\title{Simulate a magnitude production data set for both tasks and one or more
production ratios}
\usage{
simulate_gpm(ntrials, cond, param)
}
\arguments{
\item{ntrials}{A number of trials per condition.}

\item{cond}{A data frame with columns task, std, sigma, and p.
sigma: number representing the standard deviation for loudness and
brightness productions.}

\item{param}{A named vector of numbers representing the model parameters.}
}
\value{
A data frame with the simulated magnitude productions for ntrials
per condition.
}
\description{
Simulate a magnitude production data set for both tasks and one or more
production ratios
}
\examples{
cond <- data.frame(
                   task = rep(c("loud_bright", "bright_loud"), each = 6),
                   std = c(25, 34, 43, 52, 61, 70,  # loud std
                           69, 73, 77, 81, 85, 89), # bright std
                   sigma = rep(c(bright = 3,  loud = 6), each = 6),
                   p = c(1)
)
param <- create_param_set()
simulate_gpm(20, cond, param)
}
