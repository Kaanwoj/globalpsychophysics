.onLoad <- function(libname, pkgname) {
  # Check for required packages silently
  required_packages <- c("rstan", "rstantools")
  
  missing_packages <- character(0)
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      missing_packages <- c(missing_packages, pkg)
    }
  }
  
  # Store missing packages in package environment for later use
  if (length(missing_packages) > 0) {
    assign("missing_packages", missing_packages, envir = parent.env(environment()))
  }
}

.onAttach <- function(libname, pkgname) {
  # Check if there are missing packages
  if (exists("missing_packages", envir = parent.env(environment()))) {
    missing <- get("missing_packages", envir = parent.env(environment()))
    
    packageStartupMessage(
      "\n", pkgname, " requires the following packages to be installed:\n  ",
      paste(missing, collapse = ", "),
      "\n\nInstall them with:\n  install.packages(c('",
      paste(missing, collapse = "', '"), "'))",
      "\n\nFor rstan, see installation instructions at:\n  https://mc-stan.org/r-packages/\n"
    )
  } else {
    packageStartupMessage("globalpsychophysics v", utils::packageVersion(pkgname), " loaded successfully")
  }
}
