# =============================================================================
# Dependency Management for SimMSAR Project
# =============================================================================
# This file manages all package dependencies for the SimMSAR project.
# Source this file at the beginning of your scripts to ensure all required
# packages are available.
#
# Usage:
#   source("R/dependencies.R")
#
# To install missing packages:
#   source("R/dependencies.R")
#   install_missing_packages()
# =============================================================================

# Define all required packages by category
PACKAGES <- list(
  # Core MSAR functionality
  core = c(
    "NHMSAR",      # Non-homogeneous Markov-switching autoregressive models
    "netcontrol",  # Network control analysis
    "huge"         # High-dimensional undirected graph estimation (nonparanormal)
  ),

  # Data manipulation
  data = c(
    "dplyr",       # Data manipulation and transformation
    "Matrix"       # Sparse and dense matrix classes
  ),

  # Statistical modeling
  stats = c(
    "lme4",        # Linear mixed-effects models
    "lmerTest",    # Tests for mixed models
    "lmPerm",      # Permutation tests for linear models
    "dunn.test",   # Dunn's test for multiple comparisons
    "effects",     # Effect displays for linear models
    "sjPlot"       # Statistical plots and tables
  ),

  # Simulation and generation
  simulation = c(
    "graphicalVAR", # Graphical vector autoregression
    "mvtnorm",      # Multivariate normal and t distributions
    "abind"         # Combine multi-dimensional arrays
  ),

  # LASSO and regularization
  regularization = c(
    "lars",         # Least angle regression, lasso and forward stagewise
    "prettyGraphs", # Graph visualization (for LASSO)
    "glmnet",       # Lasso and elastic-net regularized GLMs (adaptive LASSO)
    "glasso"        # Graphical lasso for precision matrix estimation
  ),

  # Visualization
  visualization = c(
    "ggplot2",      # Grammar of graphics plotting
    "plotly",       # Interactive web-based graphs
    "cowplot",      # Streamlined plot theme and plot annotations
    "RColorBrewer"  # Color palettes for plots
  ),

  # Output and reporting
  output = c(
    "knitr",        # Dynamic report generation
    "kableExtra",   # Construct complex tables with kable
    "xtable"        # Export tables to LaTeX or HTML
  ),

  # Utilities
  utilities = c(
    "progress",     # Progress bars for loops
    "clue"          # Cluster ensembles, Hungarian algorithm for assignment
  )
)

# Flatten the list to get all packages
ALL_PACKAGES <- unlist(PACKAGES, use.names = FALSE)

#' Check which packages are installed
#'
#' @return Named logical vector indicating which packages are installed
check_installed_packages <- function() {
  installed <- ALL_PACKAGES %in% installed.packages()[, "Package"]
  names(installed) <- ALL_PACKAGES
  return(installed)
}

#' Install missing packages
#'
#' @param packages Character vector of package names. If NULL (default),
#'   checks all required packages.
#' @param repos Repository to install from. Default is CRAN.
#' @return Invisibly returns TRUE if all packages installed successfully
install_missing_packages <- function(packages = NULL,
                                    repos = "https://cloud.r-project.org") {
  if (is.null(packages)) {
    packages <- ALL_PACKAGES
  }

  installed <- packages %in% installed.packages()[, "Package"]
  missing <- packages[!installed]

  if (length(missing) > 0) {
    message("Installing missing packages: ", paste(missing, collapse = ", "))
    install.packages(missing, repos = repos)

    # Verify installation
    still_missing <- missing[!(missing %in% installed.packages()[, "Package"])]
    if (length(still_missing) > 0) {
      warning("Failed to install: ", paste(still_missing, collapse = ", "))
      return(invisible(FALSE))
    } else {
      message("All packages installed successfully!")
      return(invisible(TRUE))
    }
  } else {
    message("All required packages are already installed.")
    return(invisible(TRUE))
  }
}

#' Load all required packages
#'
#' @param quietly Logical. Should packages be loaded quietly? Default TRUE.
#' @return Invisibly returns a named logical vector of successfully loaded packages
load_packages <- function(quietly = TRUE) {
  message("Loading required packages...")

  loaded <- sapply(ALL_PACKAGES, function(pkg) {
    success <- suppressPackageStartupMessages(
      requireNamespace(pkg, quietly = quietly)
    )
    if (success) {
      library(pkg, character.only = TRUE, quietly = quietly)
    }
    return(success)
  })

  failed <- names(loaded)[!loaded]
  if (length(failed) > 0) {
    warning("Failed to load packages: ", paste(failed, collapse = ", "))
    message("Try installing them with: install_missing_packages()")
  } else {
    message("All packages loaded successfully!")
  }

  return(invisible(loaded))
}

#' Print package dependency summary
#'
#' @return Invisibly returns a data frame with package information
print_package_summary <- function() {
  cat("\n=== SimMSAR Package Dependencies ===\n\n")

  for (category in names(PACKAGES)) {
    cat(sprintf("%-20s (%d packages)\n",
                paste0(toupper(substring(category, 1, 1)),
                       substring(category, 2), ":"),
                length(PACKAGES[[category]])))

    for (pkg in PACKAGES[[category]]) {
      installed <- pkg %in% installed.packages()[, "Package"]
      status <- if (installed) "\u2713" else "\u2717"  # checkmark or X
      cat(sprintf("  %s %-20s", status, pkg))

      if (installed) {
        version <- as.character(packageVersion(pkg))
        cat(sprintf(" (v%s)", version))
      } else {
        cat(" [NOT INSTALLED]")
      }
      cat("\n")
    }
    cat("\n")
  }

  installed_status <- check_installed_packages()
  total <- length(installed_status)
  n_installed <- sum(installed_status)

  cat(sprintf("Total: %d/%d packages installed (%.1f%%)\n",
              n_installed, total, 100 * n_installed / total))

  invisible(data.frame(
    package = ALL_PACKAGES,
    installed = check_installed_packages(),
    stringsAsFactors = FALSE
  ))
}

# =============================================================================
# Auto-load packages when this file is sourced
# =============================================================================

# Automatically load all packages when this file is sourced
# Comment out the line below if you prefer manual control
load_packages(quietly = TRUE)

message("Dependency management loaded. Run print_package_summary() to see package status.")
