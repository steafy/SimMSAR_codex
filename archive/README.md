# Archive Directory

This directory contains files that are no longer actively used in the main workflow but are kept for reference.

## Directory Structure

### old_versions/
Deprecated version 1 files, replaced by `_2` versions in the main codebase:
- `fit.MSAR_revised.R` → replaced by `fit.MSAR_revised_2.R`
- `init.theta.MSAR_revised.R` → replaced by `init.theta.MSAR_revised_2.R`
- `as.thetaMSAR_revised.R` → replaced by `as.thetaMSAR_revised_2.R`
- `Mstep.hh.lasso.MSAR_patched.R` → replaced by `Mstep.hh.lasso.MSAR_patched_2.R`
- `Mstep.hh.reduct.MSAR_patched.R` → replaced by `Mstep.hh.reduct.MSAR_patched_2.R`

### obsolete/
Scripts that are broken or no longer functional:
- `MSARmodel_NHMSAR.R` - Old demo script that references missing files (`generate_Rseq.R`, `calc_transmat.R`) and uses deprecated version 1 files

### experimental/
Experimental or alternative implementations that were never integrated:
- `generate_Beta_2.R` - Alternative Beta generation with memoization
- `Mstep.hh.lasso.MSAR_revised.R` - Alternative M-step implementation
- `Mstep.hh.lasso_glm.MSAR_patched.R` - M-step using glmnet instead of lars
- `Mstep.hh.prune.MSAR.R` - Pruning-based M-step
- `fit.MSAR_+LASSO_alldependencies.R` - All-in-one bundled file

## Note
These files are version controlled but not part of the active codebase. They may be useful for:
- Understanding development history
- Comparing implementations
- Recovering functionality if needed

Archived: 2025-11-18
