# BIC model-selection bug in the LASSO M-step (and why parallel speedup saturated)

Investigation date: 2026-06-26. Test grid: `N={4,6,8}`, `M={1,2,3,4}`,
`density={0.25,0.5,0.75}`, `T={500,1000,2000}`, `n_ts=5` (540 fits). Machine:
6 physical / 12 logical cores, R 4.6.0, reference (single-threaded) BLAS.

## TL;DR

1. The parallel estimation's modest speedup (~2.3× on 5 workers) is **not** a
   scheduling problem — it is a hardware ceiling (memory-bandwidth bound; per-fit
   time grows ~1.4× at 6 workers, ~2.1× at 12) **plus** a pathological cost tail.
2. The cost tail was caused by a **correctness bug** in the LASSO M-step's BIC
   model-selection line (`R/estimation/mstep_hh_lasso_msar.R`), inherited verbatim
   from `NHMSAR::Mstep.hh.lasso.MSAR`.
3. **Fix applied (option 1):** correct the BIC line → **10–21× faster** on heavy
   cells **and** materially better recovery of A and Sigma.
4. **Option 2 (EBIC penalty recalibration) does not work** and was *not* applied:
   the λ-selection is structurally unable to enforce sparsity; all sparsity in the
   pipeline comes from the downstream `min_edg_val` threshold.

## The bug

`R/estimation/mstep_hh_lasso_msar.R`, inner λ-selection loop. Original line:

```r
BIC.lm[kst] = -2*sum(log(NHMSAR:::pdf.norm(
                matrix(c(wy[,id,1]), 1, length(wy)),                  # (1) wy[,id,1]  (2) length(wy)
                matrix(c(mylm[[kst]]$fitted.values), 1, length(wy)),  #               (2) length(wy)
                as.matrix(var(wy[,id,j]-mylm[[kst]]$fitted.values)))))+
              log(length(wy))*(length(w[[kst]])+2)                   # (3) length(wy)
```

`wy <- array(0, c(N.samples*(T-1), d, M))`, so for a heavy cell
(N.samples=1, T=2000, d=8, M=4) the response vector `wy[,id,j]` has length **1999**
but `length(wy)` = 1999·8·4 = **63968**.

- **(1) Wrong regime:** "observed" is `wy[,id,1]` (hardcoded regime 1), while the
  model was fit on `wy[,id,j]` and the residual variance uses `wy[,id,j]`. For any
  regime `j ≥ 2` the BIC scores regime-1 observations against regime-`j` fitted
  values — internally inconsistent in the same expression.
- **(2/3) Wrong length:** `matrix(c(v), 1, 63968)` recycles the 1999-long vector
  exactly 32× → the log-likelihood term is summed over 32× too many points
  (32× inflation), while the penalty `log()` term changes only 1.46×. The penalty
  becomes negligible relative to the inflated fit term.

### Direct evidence the selection is broken
Per-λ BIC curves from the *real* M-step (cell T=1000, D=0.5, N=6, M=2, node 1):

| nvars | 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|
| BIC (buggy), regime 2 | **48643 ← picks** | 57625 | 65347 | 66017 | 67426 |
| BIC (correct), regime 2 | 3445 | 3063 | 2852 | 2841 | **2810 ← picks** |

The buggy BIC *increases* with model size for regime 2 (because adding regime-2
structure worsens the fit to the wrong regime-1 data), so it picks the sparsest
model; the corrected BIC decreases. They select opposite models.

## The fix (option 1, APPLIED)

Use regime `j`'s own response and its actual length in all three places:

```r
n_obs = length(wy[,id,j])
BIC.lm[kst] = -2*sum(log(NHMSAR:::pdf.norm(
                matrix(c(wy[,id,j]), 1, n_obs),
                matrix(c(mylm[[kst]]$fitted.values), 1, n_obs),
                as.matrix(var(wy[,id,j]-mylm[[kst]]$fitted.values)))))+
              log(n_obs)*(length(w[[kst]])+2)
```

### Speed (A/B, identical seed/init/data, `fit_msar` MaxIter=200)
| Cell (T,D,N,M) | original | corrected | speedup | EM iters orig→fix |
|---|---|---|---|---|
| 2000,0.5,8,4 | 36.8 s | 3.15 s | 11.7× | 148 → 35 |
| 1000,0.5,8,4 | 10.2 s (failed) | 0.92 s | 11.1× | fail → 12 |
| 2000,0.75,8,3 | 21.7 s | 1.02 s | 21.2× | 200 → 6 |
| 1000,0.5,6,2 | 4.8 s | 0.49 s | 9.8× | 200 → 7 |

The 200-iteration non-converging tail was a *symptom* of the bug (over-dense
support handed to the reduced M-steps). Corrected fits converge in 6–35 iters.

### Recovery accuracy vs known-true matrices
Evaluated at the pipeline's `min_edg_val = 0.05` threshold (T∈{1000,2000}, all
N/M/density, n_ts=3, regimes matched via `match_regimes`):

| Metric | original | corrected |
|---|---|---|
| A correlation | 0.62 | **0.94** |
| A MAE | 0.186 | **0.060** |
| Sensitivity (true edges found) | 0.55 | **0.95** |
| Specificity (spurious edges avoided) | 0.844 | 0.834 |
| Edges kept (true ≈ 19) | 12.7 (under-selects) | **20.9** |
| Sigma correlation | 0.75 | **0.93** |
| Sigma Frobenius error (median) | 1.16 (116%) | **0.10 (10%)** |

On *raw* estimates the corrected version over-selects (raw specificity ~0.07), but
the spurious edges are small-weight and the 0.05 threshold prunes them, restoring
specificity to 0.83 while keeping sensitivity 0.95. The original's apparent
parsimony is just missing real edges, not good selection. `M = 1` is unaffected
(no regime mismatch).

## Option 2: EBIC penalty recalibration — TESTED, DOES NOT WORK (not applied)

Adding an Extended-BIC term `+ 2·γ·log(d)·length(w[[kst]])` to the corrected BIC:

- γ ∈ {0.5, 1, 2}: estimated support **identical** to option 1 in **0 / 468** grid
  records — no effect on any metric or edge count.
- Even **γ = 500** still selects the fully dense model (72/72 nonzero entries).

**Why:** the λ-selection refits an *unpenalized* `lm` per node and scores it by
in-sample Gaussian likelihood over ~1000–2000 (regime-weighted) points. The
likelihood gains from adding predictors dwarf any `log(n)`-scale penalty, so it
always lands on the saturated model. Consequently the "LASSO M-step" does not
itself produce sparsity — **all sparsity comes from the downstream `min_edg_val`
threshold**, in both the original and corrected code.

Genuine estimator-level sparsity would need a different mechanism, e.g.
cross-validated λ, using the LASSO/lars coefficients directly instead of refitting
an unpenalized `lm`, or stability selection (`R/estimation/bootstrap_stability.R`).

## Parallelisation note

Independent measurement: at 5 workers the measured wall (~1222 s) is ~89% of the
LPT-optimal makespan; per-replicate vs cell-bundled dispatch give identical
makespan. The workload is memory-bandwidth bound (6 physical cores), so wall-time
saturates near the physical core count — recommend `workers ≈ 6`, not more (K=8
was measured *slower* than K=6 due to hyperthread contention). The BIC fix, by
collapsing the per-fit cost, helps far more than any scheduling change.

Side note: this R 4.6.0 install segfaulted intermittently on trivial scripts
during testing — relevant for long unattended runs (a worker crash loses results).
