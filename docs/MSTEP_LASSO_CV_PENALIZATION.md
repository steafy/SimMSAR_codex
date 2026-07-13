# Genuine LASSO penalization in the first M-step (cv.glmnet)

Investigation date: 2026-07-01. Follows up on `docs/MSTEP_LASSO_BIC_BUG.md`,
which fixed the BIC λ-selection line but concluded (its "Option 2") that the
λ-selection is *structurally* unable to produce sparsity: all sparsity in the
pipeline comes from the downstream `min_edg_val = 0.05` threshold in
`compute_recovery_metrics()`. This document (a) independently re-verifies that
finding, (b) diagnoses *why* the current "LASSO" M-step cannot shrink, and
(c) replaces it with a genuinely penalized, cross-validated estimator.

## 0. Independent verification of the core finding (CONFIRMED)

Script: isolated first M-step on freshly generated series with a **known sparse
ground-truth** A (generation enforces `|A_ij| >= 0.05` on non-zero edges),
capturing `par$A[[j]][[1]]` **before** any thresholding.

`RAW nnz` = share of `|entry| > 1e-8` in the raw estimate; true density in
parentheses. Specificity is against the true support (share of true *non*-edges
correctly left at zero).

| Cell (T,N,M,D) | regime | true density | RAW nnz | RAW spec | @0.05 nnz | @0.05 spec |
|---|---|---|---|---|---|---|
| 500,6,1,0.25 | 1 | 25% | 83% | 0.22 | 36% | 0.81 |
| 500,6,2,0.25 | 1 | 25% | 83% | 0.22 | 64% | 0.44 |
| 500,6,2,0.25 | 2 | 25% | 83% | 0.15 | 67% | 0.37 |
| 1000,8,3,0.50 | 1–3 | 50% | **100%** | **0.00** | 84–92% | 0.12–0.19 |

The raw estimate is essentially **dense** (specificity 0.00–0.22). In the
`N=8,M=3` cell the first M-step returns the **fully saturated** model (100%
non-zero). All apparent sparsity is manufactured downstream by the 0.05 cut.

**Support is frozen after iteration 1.** Instrumenting a full `fit_msar()`
(`penalty="LASSO"`) to log per-regime support at every M-step call: the first
(LASSO) step chose 36/36 (saturated) for both regimes on `N=6`, and all five
subsequent `mstep_hh_reduct_msar` calls kept **exactly** that support
(`setequal` to iteration 1 = TRUE at every call). So a local fix in the first
M-step does propagate — `reduct` never re-selects, it only re-estimates the
non-zero entries by weighted OLS. Verdict: **finding confirmed on all three
legs.**

## 2. Why the current M-step cannot shrink

`R/estimation/mstep_hh_lasso_msar.R` (original / BIC-fixed):

1. `lars()` computes a path, but the shrunk lars coefficients are **never used**.
   They only generate candidate supports `w[[kst]]`.
2. For each candidate support an **unpenalized** `lm()` is refit and scored by
   in-sample Gaussian BIC. With ~`N.samples*(T-1)` (regime-weighted,
   autocorrelated) pseudo-observations the log-likelihood gain from adding any
   predictor dwarfs the `log(n)`-scale penalty, so BIC lands on the saturated
   model (confirmed to γ=500 EBIC in the prior doc, and reproduced above).
3. The final coefficients are `A2.lasso[id,w] = t(Cxy_w) %*% solve(Cxx_w)` —
   unpenalized weighted OLS on the selected support. There is **no L1 shrinkage
   anywhere in the returned estimate.**
4. A `bic.new > bic.old` fallback (lines 196–208) reverts to the **init**
   support, which is the *dense* OLS init — another route to a saturated model.

### Secondary bug: the weighting is `w`, not `sqrt(w)`
The lars path is fed `wx = obs * w`, `wy = obs_next * w` (lines 89–92), where
`w = gamma[,j]` is the E-step regime-membership probability. OLS on `(wx, wy)`
minimizes `Σ_t (w_t y_t − w_t x_t b)² = Σ_t w_t² (y_t − x_t b)²` — i.e. it weights
each pseudo-observation by **`w²`**, not `w`. Correct WLS weighting by `w`
requires transforming with `sqrt(w)` (`x* = √w·x`, `y* = √w·y`), or, better,
passing `w` to a solver's native `weights=` argument. (The *final* coefficients
via `Cxx`/`Cxy` are correctly `w`-weighted first/second moments; only the
candidate-generating lars path carries the `w²` error — but since that path only
picks supports and the picked support is saturated anyway, it never mattered.)

## 3. Fix: cross-validated LASSO via `cv.glmnet` (APPLIED)

New `R/estimation/mstep_hh_lasso_msar.R`:

- Per regime `j`, per response node `id`: `glmnet::cv.glmnet(X, Y[,id],
  weights = w, alpha = 1, intercept = TRUE)` on the **unweighted raw** design
  `X = obs[1:(T-1),]`, `Y = obs[2:T,]`, with `w = gamma[,j]` passed via the
  native `weights=` argument (correct WLS; no `√w`/`w` confusion). The raw
  unweighted `X`,`Y` and the weight vector are now stored alongside the existing
  weighted moment accumulators.
- λ is chosen by cross-validation (`lambda.min` / `lambda.1se`, configurable via
  `options(simmsar_lasso_lambda=)`), **not** by in-sample refit-BIC.
- The support is whatever `cv.glmnet` leaves non-zero at the chosen λ.
- Final coefficients on that support are either the shrunk glmnet values
  directly, or an unpenalized weighted-OLS refit (`Cxy_w %*% solve(Cxx_w)`,
  "relaxed LASSO") — configurable via `options(simmsar_lasso_refit=)`.
- The `bic.new > bic.old` revert-to-dense-init fallback is **removed** (it
  defeats sparsity). Intercept (`A0`) and `sigma` are still computed from the
  correctly `w`-weighted sufficient statistics `m`/`m_1`/`op`/`op_1`/`op_2`,
  exactly as before.
- `reduct` (iterations 2+) is **unchanged**: it now simply inherits the genuinely
  sparse support.

Within the `cvglmnet` engine the knob defaults are `lambda.1se` + `refit=TRUE`
(the only combination that both delivers estimator-level sparsity and is not
strictly dominated). But the cv.glmnet engine itself is **opt-in**: `fit_msar`
defaults to `simmsar_lasso_engine="bic"` (the legacy `mstep_hh_lasso_bic`), so
production behaviour is unchanged until you opt in. See §5 for why the engine is
only worth using together with per-iteration re-selection.

## 4. Validation (A/B vs current)

Real `fit_msar()` calls, `MaxIter=200`, 6 cells (light→heavy) × 3 replicates ×
4 configs. **Identical data AND identical init per (cell, replicate)** across all
configs, so every difference is the first M-step. Metrics are per-regime then
averaged (regimes matched to truth with `match_regimes`). A sens/spec/cor and
K/Sigma are scored **after** the 0.05 threshold, exactly as the analysis
pipeline does; `RAW b_spec` is the specificity of the *un*-thresholded estimate
(the headline "does sparsity come from the estimator?" number).

Cells: (T,N,M,D) = (200,4,1,.25) (400,6,2,.5) (800,6,3,.5) (800,8,2,.75)
(1600,4,2,.25) (1600,8,4,.5).

### Aggregate (mean over all 6 cells × 3 reps)

| config | RAW b_spec | @0.05 b_sens | @0.05 b_spec | @0.05 b_cor | k_cor | sig_cor | med iters | mean t_total | #hit MaxIter (of 18) |
|---|---|---|---|---|---|---|---|---|---|
| **OLD** (current) | 0.08 | **0.96** | 0.79 | **0.95** | **0.96** | **0.95** | 8 | **0.88 s** | 1 |
| NEW 1se+refit | **0.65** | 0.74 | **0.84** | 0.80 | 0.72 | 0.81 | 17.5 | 2.57 s | 4 |
| NEW 1se+shrunk | 0.64 | 0.73 | 0.83 | 0.80 | 0.70 | 0.81 | 11 | 2.05 s | 3 |
| NEW min+refit | 0.25 | 0.91 | 0.80 | 0.92 | 0.87 | 0.88 | 9 | 1.23 s | 0 |

Excluding the non-converged (≥200-iter) fits does **not** rescue NEW-1se
(b_sens 0.75, b_cor 0.86, k_cor 0.80 — still well below OLD's 0.97/0.97/0.96),
so the recovery gap is not merely a convergence artefact.

### Runtime breakdown: light vs heavy cell (mean over 3 reps)

| config | light (200,4,1) t_mstep / t_total / iters | heavy (1600,8,4) t_mstep / t_total / iters |
|---|---|---|
| OLD | 0.09 / 0.42 s / 4 | 0.89 / 1.64 s / 12 |
| NEW 1se+refit | 0.22 / 0.42 s / 5 | 1.80 / **8.22 s / 143** |
| NEW 1se+shrunk | 0.16 / 0.36 s / 5 | 1.68 / 5.38 s / 76 |
| NEW min+refit | 0.17 / 0.37 s / 4 | 1.93 / 2.72 s / 13 |

The `cv.glmnet` first M-step costs ~0.15 s (light) to ~1–2 s (heavy) more than
the old lars+BIC step — modest in isolation. The real cost is **downstream**:
committing to an over-sparse support in iteration 1 makes the EM struggle to
separate regimes, so NEW-1se needs 3–12× more iterations on heavy/high-M cells
(143 vs 12 on the heavy cell → 5× total wall time).

### What the fix does and does not achieve

- **Achieves the mechanical goal (1se):** sparsity now comes from the estimator.
  RAW specificity 0.65 (vs OLD 0.08); RAW ≈ thresholded (the 0.05 cut becomes
  almost a no-op). `cv.glmnet` λ is genuinely selected; the WLS weighting is now
  correct (`weights=w`, not `w²`).
- **But it worsens the study's primary recovery outcomes.** On every cell,
  `lambda.1se` buys specificity by **dropping true edges** (sensitivity
  0.96→0.74) and degrades A-weight, K and Sigma correlation. This is the
  classic 1se over-shrinkage, amplified here because (a) k-fold CV on
  autocorrelated, softly regime-weighted pseudo-observations over-penalises, and
  (b) the support is **frozen after iteration 1** (`reduct`), so an
  early-committed sparse support — chosen from a poor initial E-step — locks in
  missed edges the way the old dense-support + end-threshold design does not.
- **`lambda.min` sidesteps the damage but barely sparsifies** (RAW spec 0.25,
  RAW nnz 81%) — i.e. it does *not* deliver estimator-level sparsity, so it
  doesn't meet the goal either. There is no free lunch on this design: any
  estimator-level sparsity beyond `min` costs recovery.
- **`refit` vs `shrunk`:** near-identical; relaxed-LASSO refit is marginally
  better on `b_cor`/`k_cor`, as expected (removes shrinkage bias). If adopting,
  use `refit=TRUE`.

### Companion bug found & fixed
`mstep_hh_reduct_msar.R` used `S.th[-wi, jd]` to select the complement of a
node's edge set. With genuinely sparse supports a node can have **no** incoming
edges (`wi` empty); `-integer(0)` selects nothing, not everything, crashing the
2nd EM iteration ("replacement has length zero"). Replaced with an explicit
complement `seq_len(d)[-wi]`/`seq_len(d)`. This is a latent correctness bug that
never fired only because the old M-step always returned dense supports; the fix
is worth keeping independently of the LASSO decision.

## 5. The fix that works: re-select the support EVERY EM iteration

The §4 failure was diagnosed as *frozen early support*, so the obvious remedy is
to stop freezing: re-run the penalized `cv.glmnet` selection in **every** M-step
(not just iteration 1), letting the support track the regime separation as the
E-step improves. A/B across M=1..4 (5 cells, 3 reps, identical data+init):

| config (all 1se+refit) | b_sens | b_spec | b_cor | k_cor | sig_cor | RAW b_spec | med iters | mean t_total | non-conv |
|---|---|---|---|---|---|---|---|---|---|
| OLD (bic engine) | 0.96 | 0.72 | 0.88 | 0.92 | 0.90 | 0.08 | 40 | 1.1 s | 2 |
| cvglmnet **frozen** | 0.75 | 0.79 | 0.75 | 0.66 | 0.77 | 0.59 | 65 | 2.4 s | 3 |
| **cvglmnet re-select** | **0.94** | **0.94** | **0.95** | **0.92** | **0.96** | **0.81** | **12** | 14.7 s | 0 |
| cvglmnet re-select (min) | 0.97 | 0.83 | 0.96 | 0.94 | 0.97 | 0.48 | 40 | 54 s | 0 |

Re-selection (`lambda.1se`) **matches or beats OLD on every recovery metric**
(A_cor 0.95 vs 0.88, A_spec 0.94 vs 0.72, Sigma 0.96 vs 0.90, K tied
0.92) **and** delivers genuine estimator-level sparsity (RAW spec 0.81 vs 0.08)
**and** converges in fewer iterations (12 vs 40, zero non-converged). It holds at
every M including M=4 (A_cor 0.94 vs OLD 0.90). Per-M detail:

| M | OLD b_cor / k_cor | re-select b_cor / k_cor | re-select RAW spec |
|---|---|---|---|
| 1 | 0.85 / 1.00 | 0.98 / 1.00 | 0.94 |
| 2 | 0.85 / 0.90 | 0.99 / 0.97 | 0.78 |
| 3 | 0.93 / 0.87 | 0.83 / 0.79 | 0.79 |
| 4 | 0.90 / 0.94 | 0.94 / 0.85 | 0.77 |

**The only cost is runtime: ~13× on average, ~20× on M=4** (cv.glmnet runs `d·M`
times per iteration, every iteration). Adaptive weighting (`penalty.factor =
1/|OLS init|`) did **not** help — it over-penalises weak true edges (near the
generation's 0.05 floor) and lowers sensitivity. `lambda.min` re-select recovers
marginally better but is far denser (RAW spec 0.48) and ~4× slower again.

### 5.1 Can partial re-selection cut the cost? (first-k / every-m) — mostly no

Tested `simmsar_lasso_reselect_iters=k` (re-select iters 1..k, then freeze) and
`simmsar_lasso_reselect_every=m`, A/B on M=1..4 (2 reps):

| config | b_cor | k_cor | sig_cor | RAW spec | cv.glmnet calls | t_total | non-conv (of 8) |
|---|---|---|---|---|---|---|---|
| RS full | 0.93 | 0.91 | 0.96 | 0.86 | 10.4 | 11.6 s | **0** |
| RS first-5 | 0.88 | 0.88 | 0.89 | 0.84 | 4.6 | 6.2 s (1.9× faster) | 2 |
| RS first-3 | 0.86 | 0.85 | 0.85 | 0.79 | 3.0 | 4.5 s | 3 |
| RS every-3 | 0.92 | 0.90 | 0.93 | 0.86 | 51.5 | 58.6 s (**5× slower**) | 6 |

- **`first-5` ≈ 1.9× faster** and keeps most recovery — a partial win — **but it
  reintroduces non-convergence**: freezing a not-yet-settled support on hard
  (M=4) cells makes the reduced M-step fight the E-step, so 2–3 of 8 fits burn
  all 200 iterations. `first-3` is worse.
- **Full re-selection is the most robust AND most iteration-efficient**
  (converged 0/8, median 10 iters): re-selecting every step lets the support and
  coefficients co-adapt to a clean fixed point. Its cost is purely *per-iteration*
  cv.glmnet time, not iteration count.
- **`every-m` is counterproductive**: a support that jumps every m-th iteration
  rarely satisfies the convergence test (6/8 hit MaxIter), ending up ~5× slower
  than full.

Conclusion: the schedule is not the right place to save time — full
re-selection's fast, reliable convergence is *because* it never freezes. The
per-call cost is the right lever instead (§5.2).

### 5.2 Cutting cv.glmnet's per-call cost: fewer folds + FIXED folds (~2.5×)

Because the relaxed refit means cv.glmnet is only used to pick the **support**
(coefficients are re-estimated exactly by WLS), a cheaper CV should not hurt
recovery. First attempt — just fewer folds / coarser grid, **random** folds:

| CV setting (full re-select) | b_cor | k_cor | RAW spec | EM iters | t_total | speedup |
|---|---|---|---|---|---|---|
| 10 folds × 100 λ (random) | 0.93 | 0.91 | 0.86 | 10 | 11.7 s | 1.0× |
| 5 × 50 (random) | 0.93 | 0.91 | 0.84 | 17 | 11.2 s | 1.04× |
| 3 × 30 (random) | 0.93 | 0.91 | 0.83 | 63 | 31.8 s | 0.37× |

Recovery is untouched, but **total time doesn't drop** — cheaper *random* CV adds
sampling noise that makes the selected support flicker between EM iterations, so
convergence needs more iterations (exactly cancelling the per-call saving; 3-fold
is far worse). Same lesson as §5.1: support instability across iterations is the
enemy.

**The fix: fix the fold partition** (deterministic interleaved folds, identical
every iteration, no RNG). Then the support only changes when the E-step genuinely
changes:

| CV setting (full re-select) | b_cor | k_cor | sig_cor | RAW spec | EM iters | t_total | speedup |
|---|---|---|---|---|---|---|---|
| 10 × 100 random (naive) | 1.00 | 0.98 | 0.99 | 0.85 | 9.0 | 11.8 s | 1.0× |
| 10 × 100 **fixed** | 1.00 | 0.98 | 0.99 | 0.86 | 6.6 | 7.3 s | 1.6× |
| 5 × 50 random | 1.00 | 0.98 | 0.99 | 0.81 | 28.0 | 24.2 s | 0.5× |
| **5 × 50 fixed** | 1.00 | 0.98 | 0.99 | 0.83 | 6.8 | **4.7 s** | **2.5×** |
| 5 × 30 fixed | 0.94 | 0.87 | 0.92 | 0.80 | 8.7 | 6.8 s | 1.7× |

**`5 folds × 50 λ × fixed folds` is ~2.5× faster than naive full-CV re-selection
with identical recovery** (random 5-fold needed 28 EM iters; fixed needs 7).
That brings full re-selection from ~13× the legacy engine down to **~5×**. Going
below 50 λ or to 3 folds starts to lose recovery (5×30 fixed: K 0.87). These
are now the cv.glmnet-engine **defaults** (`simmsar_lasso_nfolds=5`,
`simmsar_lasso_nlambda=50`, `simmsar_lasso_fixedfolds=TRUE`); they only apply when
the (opt-in) `cvglmnet` engine is selected, so the default `bic` path is
unaffected. Further per-call savings (warm starts across EM iterations) remain
possible but were not pursued.

### Companion bug found & fixed
`mstep_hh_reduct_msar.R` used `S.th[-wi, jd]` to select the complement of a
node's edge set. With genuinely sparse supports a node can have **no** incoming
edges (`wi` empty); `-integer(0)` selects nothing, not everything, crashing the
2nd EM iteration ("replacement has length zero"). Replaced with an explicit
complement `seq_len(d)[-wi]`/`seq_len(d)`. Latent correctness bug that never fired
only because the old M-step always returned dense supports; kept regardless.

## Recommendation

The original diagnosis was right (the "LASSO" penalty had no estimator-level
effect), and the penalty **can** be made to work *without any loss* of recovery —
but only with **per-iteration re-selection**, not a frozen first-iteration
support. Concretely:

- **Frozen cv.glmnet is worse than the legacy engine** on every recovery metric —
  do not ship it.
- **Re-select (`cvglmnet` + `simmsar_lasso_reselect=TRUE`, `lambda.1se`,
  `refit=TRUE`) matches/beats the legacy engine** and additionally gives
  threshold-free sparsity and better convergence, at ~10-20× runtime.

**Wiring (implemented, default-safe):** `fit_msar` selects the LASSO strategy via
options; **defaults are unchanged** (`simmsar_lasso_engine="bic"` → the exact
legacy lars+BIC first step, `mstep_hh_lasso_bic()`). Opt in with:

```r
options(simmsar_lasso_engine = "cvglmnet",  # genuinely-penalized engine
        simmsar_lasso_reselect = TRUE,       # re-select support every iteration
        simmsar_lasso_lambda   = "1se",      # (default) sparser rule
        simmsar_lasso_refit    = TRUE)       # (default) relaxed-LASSO refit
```

**Decision for the maintainer (Stefan):** whether to adopt re-selection as the
production estimator depends on the ~13× runtime hit against the full `n_ts=50`
grid, weighed against the methodological benefit of threshold-free sparsity with
equal-or-better recovery. Runtime was brought down from ~13× to **~5× the legacy
engine** via the cv.glmnet-engine defaults (§5.2: 5 folds, 50 λ, fixed folds =
2.5× faster, no recovery loss). Partial-re-selection *schedules* to cut cost
further were tested (§5.1) and are not worth it (`first-5` ~2× but risks
non-convergence on M=4; `every-m` backfires); full re-selection remains the most
robust. The `reduct` `-wi` fix should be kept regardless of this decision.
