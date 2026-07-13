# =============================================================================
# Dependency management for SimMSAR
# =============================================================================
# Source this file at the start of a session. It CHECKS that the required
# packages are installed and, if any are missing, prints the exact install
# command(s) and stops -- it never installs anything automatically.
#
#   source("R/dependencies.R")     # checks + loads, or stops with instructions
#   print_package_summary()        # optional: per-package installed/version table
# =============================================================================

# Archived CRAN package: not installable via a plain install.packages(). NHMSAR
# provides the MSAR EM machinery adapted under R/estimation/.
ARCHIVED_PACKAGES <- c(
  NHMSAR = "https://cran.r-project.org/src/contrib/Archive/NHMSAR/NHMSAR_1.19.tar.gz"
)

# CRAN packages actually used by the code (verified by usage, not guessed).
PACKAGES <- list(
  estimation  = c("glmnet", "lars"),                       # LASSO / lars M-steps
  generation  = c("mvtnorm", "abind"),                     # residual draws, arrays
  data        = c("dplyr", "tibble", "stringr"),           # wrangling
  parallel    = c("future", "future.apply", "progressr"),  # parallel backend + progress
  models      = c("lme4", "lmerTest", "emmeans",           # mixed models + contrasts
                  "performance", "dunn.test"),
  viz         = c("ggplot2", "cowplot", "viridisLite",     # plots
                  "systemfonts"),
  tables      = c("kableExtra", "xtable"),                 # LaTeX/HTML tables
  utils       = c("clue")                                  # Hungarian assignment
)

# Optional: enables rounded facet-card styling in the plots (GitHub-only). The
# plotting code falls back gracefully if it is absent, so it is not required.
#   remotes::install_github("teunbrand/elementalist")
OPTIONAL_PACKAGES <- c("elementalist")

ALL_PACKAGES <- c(names(ARCHIVED_PACKAGES), unlist(PACKAGES, use.names = FALSE))

# Logical vector: is each required package installed?
check_installed_packages <- function() {
  installed <- ALL_PACKAGES %in% rownames(installed.packages())
  names(installed) <- ALL_PACKAGES
  installed
}

# Check dependencies and STOP (with actionable instructions) if any are missing.
# Never installs anything. Called automatically when this file is sourced.
check_dependencies <- function() {
  have <- rownames(installed.packages())

  missing_cran     <- setdiff(unlist(PACKAGES, use.names = FALSE), have)
  missing_archived <- setdiff(names(ARCHIVED_PACKAGES), have)

  if (length(missing_cran) == 0 && length(missing_archived) == 0) {
    return(invisible(TRUE))
  }

  msg <- c("Missing required packages -- install them and re-source, then retry.\n")

  if (length(missing_cran) > 0) {
    msg <- c(msg, "\nCRAN packages:\n",
             sprintf('  install.packages(c(%s))\n',
                     paste(sprintf('"%s"', missing_cran), collapse = ", ")))
  }

  if (length(missing_archived) > 0) {
    urls <- ARCHIVED_PACKAGES[missing_archived]
    msg <- c(msg,
             "\nArchived CRAN package(s) -- NOT on the current CRAN, install from the archive:\n",
             "  # install.packages('remotes')  # if needed\n",
             paste0(sprintf('  remotes::install_url("%s")\n', urls), collapse = ""))
  }

  stop(paste0(msg, collapse = ""), call. = FALSE)
}

# Attach all required packages (assumes check_dependencies() has passed).
load_packages <- function(quietly = TRUE) {
  invisible(lapply(ALL_PACKAGES, function(pkg) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE, quietly = quietly))
  }))
}

# Print a per-package installed/version table (optional diagnostic).
print_package_summary <- function() {
  have <- rownames(installed.packages())
  report_block <- function(title, pkgs) {
    cat(sprintf("%-20s (%d packages)\n", title, length(pkgs)))
    for (pkg in pkgs) {
      installed <- pkg %in% have
      status <- if (installed) "✓" else "✗"
      cat(sprintf("  %s %-18s", status, pkg))
      cat(if (installed) sprintf(" (v%s)\n", as.character(packageVersion(pkg))) else " [NOT INSTALLED]\n")
    }
    cat("\n")
  }
  cat("\n=== SimMSAR package dependencies ===\n\n")
  report_block("Archived (CRAN):", names(ARCHIVED_PACKAGES))
  for (category in names(PACKAGES)) {
    report_block(paste0(toupper(substring(category, 1, 1)), substring(category, 2), ":"),
                 PACKAGES[[category]])
  }
  report_block("Optional:", OPTIONAL_PACKAGES)
  n_ok <- sum(check_installed_packages())
  cat(sprintf("Total required: %d/%d installed.\n", n_ok, length(ALL_PACKAGES)))
  invisible(data.frame(package = ALL_PACKAGES,
                       installed = check_installed_packages(),
                       row.names = NULL))
}

# On source: verify dependencies (stops with instructions if missing), then load.
check_dependencies()
load_packages(quietly = TRUE)
