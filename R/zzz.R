# Copyright (C) 2015, 2016, 2017 Trustees of Columbia University

.onLoad <- function(libname, pkgname) { # nocov start
  modules <- paste0("stan_fit4", names(stanmodels), "_mod")
  for (m in modules) loadModule(m, what = TRUE)
} # nocov end

.onAttach <- function(...) {
  stansurvLib <- dirname(system.file(package = "stansurv"))
  pkgdesc <- suppressWarnings(utils::packageDescription("stansurv", lib.loc = stansurvLib))
  if (length(pkgdesc) > 1) {
    builddate <- gsub(';.*$', '', pkgdesc$Packaged)
    packageStartupMessage(paste("stansurv (Version ", pkgdesc$Version, ", packaged: ", builddate, ")", sep = ""))
  }
  packageStartupMessage("- Do not expect the default priors to remain the same in future rstanarm versions.")
  packageStartupMessage("Thus, R scripts should specify priors explicitly, even if they are just the defaults.")
  packageStartupMessage("- For execution on a local, multicore CPU with excess RAM we recommend calling")
  packageStartupMessage("options(mc.cores = parallel::detectCores())")
  packageStartupMessage("- Plotting theme set to bayesplot::theme_default().")
  ggplot2::theme_set(bayesplot::theme_default())
}

