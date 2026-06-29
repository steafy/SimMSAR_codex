# =============================================================================
# Centralised output paths for the analysis pipeline
# =============================================================================
# All analysis artefacts are written under output/ in two subfolders:
#   * results : tables / text (Results_*_mixed.{html,tex}, comparison_all.tex,
#               descriptive_statistics_table.{html,tex}, feasibility_model.tex)
#   * plots   : figures (line / sens-spec / coefficient plots, the diagnostics
#               PDF) -- everything produced by ggsave()/pdf().
# Use results_file()/plots_file() at every write site instead of a bare
# filename so the destinations stay in one place and the directories are
# created on demand. Paths are relative to the project root (the pipeline is
# always run with cwd = project root).

RESULTS_DIR <- file.path("output", "results")
PLOTS_DIR   <- file.path("output", "plots")

# Build a path under output/results, creating the directory if needed.
results_file <- function(...) {
  if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)
  file.path(RESULTS_DIR, ...)
}

# Build a path under output/plots, creating the directory if needed.
plots_file <- function(...) {
  if (!dir.exists(PLOTS_DIR)) dir.create(PLOTS_DIR, recursive = TRUE)
  file.path(PLOTS_DIR, ...)
}
