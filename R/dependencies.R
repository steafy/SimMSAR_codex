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

# Archived CRAN packages (no longer available via install.packages).
# Installed from the CRAN archive .tar.gz via devtools::install_url().
ARCHIVED_PACKAGES <- c(
  NHMSAR     = "https://cran.r-project.org/src/contrib/Archive/NHMSAR/NHMSAR_1.19.tar.gz",
  netcontrol = "https://cran.r-project.org/src/contrib/Archive/netcontrol/netcontrol_0.1.tar.gz"
)

# Define all required packages by category (CRAN-available only)
PACKAGES <- list(
  # Core MSAR functionality
  core = c(
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
    "lmtest",      # Diagnostic tests for linear models
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
    "clue",         # Cluster ensembles, Hungarian algorithm for assignment
    "future",       # Unified parallel-execution backend (multisession etc.)
    "future.apply", # apply-family functions (future_lapply) over future plans
    "progressr"     # Cross-process progress reporting, future-aware
  )
)

# Flatten the list to get all packages (CRAN + archived)
ALL_PACKAGES <- c(names(ARCHIVED_PACKAGES), unlist(PACKAGES, use.names = FALSE))

check_installed_packages <- function() {
  installed <- ALL_PACKAGES %in% installed.packages()[, "Package"]
  names(installed) <- ALL_PACKAGES
  return(installed)
}

install_missing_packages <- function(packages = NULL,
                                    repos = "https://cloud.r-project.org") {
  if (is.null(packages)) {
    packages <- ALL_PACKAGES
  }

  installed <- packages %in% installed.packages()[, "Package"]
  missing <- packages[!installed]

  if (length(missing) == 0) {
    message("All required packages are already installed.")
    return(invisible(TRUE))
  }

  # Split into archived vs CRAN packages
  missing_archived <- missing[missing %in% names(ARCHIVED_PACKAGES)]
  missing_cran     <- missing[!missing %in% names(ARCHIVED_PACKAGES)]

  # Install CRAN packages
  if (length(missing_cran) > 0) {
    message("Installing CRAN packages: ", paste(missing_cran, collapse = ", "))
    install.packages(missing_cran, repos = repos)
  }

  # Install archived packages from CRAN archive .tar.gz via devtools
  if (length(missing_archived) > 0) {
    if (!requireNamespace("devtools", quietly = TRUE)) {
      message("Installing devtools (needed to rebuild archived packages)...")
      install.packages("devtools", repos = repos)
    }
    for (pkg in missing_archived) {
      url <- ARCHIVED_PACKAGES[[pkg]]
      message("Installing archived package '", pkg, "' from ", url)
      devtools::install_url(url)
    }
  }

  # Verify installation
  still_missing <- missing[!(missing %in% installed.packages()[, "Package"])]
  if (length(still_missing) > 0) {
    warning("Failed to install: ", paste(still_missing, collapse = ", "))
    return(invisible(FALSE))
  } else {
    message("All packages installed successfully!")
    return(invisible(TRUE))
  }
}

load_packages <- function(quietly = TRUE) {
  message("Loading required packages...")

  installed_pkgs <- installed.packages()[, "Package"]

  loaded <- sapply(ALL_PACKAGES, function(pkg) {
    if (!pkg %in% installed_pkgs) return(FALSE)
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
    message("Not yet installed: ", paste(failed, collapse = ", "))
    message("Run install_missing_packages() to install them.")
  } else {
    message("All packages loaded successfully!")
  }

  return(invisible(loaded))
}

print_package_summary <- function() {
  cat("\n=== SimMSAR Package Dependencies ===\n\n")

  # Show archived packages first
  if (length(ARCHIVED_PACKAGES) > 0) {
    cat(sprintf("%-20s (%d packages)\n", "Archived (CRAN):", length(ARCHIVED_PACKAGES)))
    for (pkg in names(ARCHIVED_PACKAGES)) {
      installed <- pkg %in% installed.packages()[, "Package"]
      status <- if (installed) "\u2713" else "\u2717"
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
# Auto-install and load packages when this file is sourced
# =============================================================================

install_missing_packages()
load_packages(quietly = TRUE)
