# =============================================================================
# migrate_rds_names.R
# =============================================================================
# One-off migration for pre-cleanup .rds files (MSAR_models_*.rds AND
# Timeseries_data_*.rds).
#
# The publication cleanup renamed the internal network notation from Beta/Kappa
# to A/K/AC (temporal network A, contemporaneous/precision network K, average
# controllability AC). Datasets computed BEFORE that rename still carry the old
# names. This script renames them so an existing .rds loads/estimates cleanly
# under the new code, WITHOUT recomputing anything: only names change, values
# and structure are untouched. The original file is never modified; a copy is
# written next to it as <name>_renamed.rds.
#
# Two schemas are handled (auto-detected; a file is usually one or the other):
#
# 1. MSAR_models_*.rds (estimation output) -- rename data-frame COLUMNS:
#      orig_Beta -> orig_A      est_Beta -> est_A
#      orig_Kappa -> orig_K     est_Kappa -> est_K
#      orig_Beta_ac -> orig_AC  est_Beta_ac -> est_AC
#      zero_var_beta -> zero_var_A   partial_zero_var_beta -> partial_zero_var_A
#
# 2. Timeseries_data_*.rds (generation output) -- rename the ELEMENT names inside
#    the `regime_dynamics` list-column (each regime's dynamics list):
#      Beta -> A, Beta_sd -> A_sd, Beta_strength -> A_strength,
#      Beta_mean_strength -> A_mean_strength, Beta_density -> A_density,
#      Beta_weighted_density -> A_weighted_density, Beta_stability -> A_stability,
#      Beta_ac -> AC, kappa -> K, kappa_pos.definit -> K_pos.definit,
#      kappa_thresh -> K_thresh, kappa_density -> K_density,
#      kappa_weighted_density -> K_weighted_density
#
# NOT renamed: Sigma/sigma (residual covariance) and mu; Cohen's-kappa fields
# (regime-sequence recovery, RQ4) in the MSAR "fit_results" attribute.
#
# NOTE: CONFIG_*.rds contains no Beta/Kappa names and needs no migration.
#
# Usage:
#   Rscript scripts/migrate_rds_names.R path/to/OLD.rds
# or set INPUT_PATH below and source() the file.
# =============================================================================

# --- Input path: command-line argument, or edit this fallback ----------------
INPUT_PATH <- "path/to/OLD.rds"

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) INPUT_PATH <- args[[1]]

if (!file.exists(INPUT_PATH)) {
  stop("Input file not found: ", INPUT_PATH,
       "\nPass the path as an argument or edit INPUT_PATH at the top of the script.")
}

# --- Exact old -> new name maps ----------------------------------------------
# MSAR_models column names:
column_map <- c(
  orig_Beta             = "orig_A",
  est_Beta              = "est_A",
  orig_Kappa            = "orig_K",
  est_Kappa             = "est_K",
  orig_Beta_ac          = "orig_AC",
  est_Beta_ac           = "est_AC",
  zero_var_beta         = "zero_var_A",
  partial_zero_var_beta = "partial_zero_var_A"
)

# Timeseries regime_dynamics element names:
regime_map <- c(
  Beta                  = "A",
  Beta_sd               = "A_sd",
  Beta_strength         = "A_strength",
  Beta_mean_strength    = "A_mean_strength",
  Beta_density          = "A_density",
  Beta_weighted_density = "A_weighted_density",
  Beta_stability        = "A_stability",
  Beta_ac               = "AC",
  kappa                 = "K",
  kappa_pos.definit     = "K_pos.definit",
  kappa_thresh          = "K_thresh",
  kappa_density         = "K_density",
  kappa_weighted_density = "K_weighted_density"
)

# Rename only exact matches present in `nms`; leave everything else untouched.
rename_names <- function(nms, map) {
  hit <- nms %in% names(map)
  nms[hit] <- map[nms[hit]]
  nms
}

# --- Load --------------------------------------------------------------------
obj <- readRDS(INPUT_PATH)

if (!is.data.frame(obj)) {
  stop("Expected a data.frame/tibble in ", INPUT_PATH,
       " but got an object of class: ", paste(class(obj), collapse = ", "))
}

changed <- character(0)

# --- (1) MSAR_models: rename data-frame columns ------------------------------
old_names <- names(obj)
names(obj) <- rename_names(old_names, column_map)
renamed_cols <- old_names[old_names != names(obj)]
if (length(renamed_cols) > 0) {
  changed <- c(changed, sprintf("column %s -> %s",
                                renamed_cols, column_map[renamed_cols]))
}

# Also rename inside the "fit_results" attribute if it carries any mapped column
# (in practice it does not; kept for self-consistency).
fr <- attr(obj, "fit_results")
if (is.data.frame(fr)) {
  names(fr) <- rename_names(names(fr), column_map)
  attr(obj, "fit_results") <- fr
}

# --- (2) Timeseries_data: rename nested regime_dynamics element names ---------
if ("regime_dynamics" %in% names(obj)) {
  n_touched <- 0L
  obj$regime_dynamics <- lapply(obj$regime_dynamics, function(fit) {
    if (!is.list(fit)) return(fit)
    lapply(fit, function(reg) {           # reg = one regime's dynamics list
      if (is.list(reg) && !is.null(names(reg))) {
        new <- rename_names(names(reg), regime_map)
        if (!identical(new, names(reg))) n_touched <<- n_touched + 1L
        names(reg) <- new
      }
      reg
    })
  })
  if (n_touched > 0) {
    changed <- c(changed, sprintf(
      "regime_dynamics element names (%d regime lists): %s",
      n_touched,
      paste(sprintf("%s->%s", names(regime_map), regime_map), collapse = ", ")))
  }
}

# --- Report ------------------------------------------------------------------
if (length(changed) == 0) {
  message("No old Beta/Kappa names found; nothing to rename. ",
          "The file may already use the A/K/AC notation.")
} else {
  message("Renamed:\n  ", paste(changed, collapse = "\n  "))
}

# --- Write migrated copy (never overwrite the original) ----------------------
out_path <- sub("\\.rds$", "_renamed.rds", INPUT_PATH, ignore.case = TRUE)
if (identical(out_path, INPUT_PATH)) out_path <- paste0(INPUT_PATH, "_renamed.rds")
saveRDS(obj, out_path)
message("Wrote migrated dataset to: ", out_path)
message("Original file left unchanged: ", INPUT_PATH)
