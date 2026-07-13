# =============================================================================
# migrate_rds_names.R
# =============================================================================
# One-off migration for pre-cleanup MSAR_models_*.rds files.
#
# The publication cleanup renamed the internal network notation from Beta/Kappa
# to A/K/AC (temporal network A, contemporaneous/precision network K, average
# controllability AC). Datasets computed BEFORE that rename still carry the old
# list-column names. This script renames those columns in place so an existing
# .rds loads cleanly under the new code, WITHOUT recomputing anything: only
# column names change, values and structure are untouched.
#
# Sigma columns (orig_Sigma / est_Sigma) are the residual covariance and are
# NOT renamed. Cohen's-kappa columns (regime-sequence recovery, RQ4) live in
# the analysis-time "fit_results" attribute under names like cohens_kappa and
# are likewise left alone.
#
# Column mapping (old -> new):
#   orig_Beta            -> orig_A
#   est_Beta             -> est_A
#   orig_Kappa           -> orig_K
#   est_Kappa            -> est_K
#   orig_Beta_ac         -> orig_AC
#   est_Beta_ac          -> est_AC
#   zero_var_beta        -> zero_var_A
#   partial_zero_var_beta-> partial_zero_var_A
#
# Usage:
#   Rscript scripts/migrate_rds_names.R path/to/MSAR_models_OLD.rds
# or set INPUT_PATH below and source() the file. The migrated copy is written
# next to the input as <name>_renamed.rds; the original is never modified.
# =============================================================================

# --- Input path: command-line argument, or edit this fallback ----------------
INPUT_PATH <- "path/to/MSAR_models_OLD.rds"

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) INPUT_PATH <- args[[1]]

if (!file.exists(INPUT_PATH)) {
  stop("Input file not found: ", INPUT_PATH,
       "\nPass the path as an argument or edit INPUT_PATH at the top of the script.")
}

# --- Exact old -> new column-name mapping ------------------------------------
name_map <- c(
  orig_Beta             = "orig_A",
  est_Beta              = "est_A",
  orig_Kappa            = "orig_K",
  est_Kappa             = "est_K",
  orig_Beta_ac          = "orig_AC",
  est_Beta_ac           = "est_AC",
  zero_var_beta         = "zero_var_A",
  partial_zero_var_beta = "partial_zero_var_A"
)

# Rename only exact matches present in `nms`; leave everything else untouched.
rename_exact <- function(nms) {
  hit <- nms %in% names(name_map)
  nms[hit] <- name_map[nms[hit]]
  nms
}

# --- Load, rename, report ----------------------------------------------------
obj <- readRDS(INPUT_PATH)

if (!is.data.frame(obj)) {
  stop("Expected a data.frame/tibble in ", INPUT_PATH,
       " but got an object of class: ", paste(class(obj), collapse = ", "))
}

old_names <- names(obj)
names(obj) <- rename_exact(old_names)
renamed_cols <- old_names[old_names != names(obj)]

# Also rename inside the "fit_results" attribute if it happens to carry any of
# the mapped columns (in practice it does not, but this keeps the migration
# self-consistent should an older schema have stored them there).
fr <- attr(obj, "fit_results")
if (is.data.frame(fr)) {
  names(fr) <- rename_exact(names(fr))
  attr(obj, "fit_results") <- fr
}

if (length(renamed_cols) == 0) {
  message("No old Beta/Kappa column names found; nothing to rename. ",
          "The file may already use the A/K/AC notation.")
} else {
  message("Renamed columns: ",
          paste(sprintf("%s -> %s", renamed_cols, name_map[renamed_cols]),
                collapse = ", "))
}

# --- Write migrated copy (never overwrite the original) ----------------------
out_path <- sub("\\.rds$", "_renamed.rds", INPUT_PATH, ignore.case = TRUE)
if (identical(out_path, INPUT_PATH)) out_path <- paste0(INPUT_PATH, "_renamed.rds")
saveRDS(obj, out_path)
message("Wrote migrated dataset to: ", out_path)
message("Original file left unchanged: ", INPUT_PATH)
