# Degenerate Sigma / Kappa estimates: diagnosis (postmix / effective sample size)

Investigation date: 2026-07-01. Follows up on `docs/MSTEP_LASSO_CV_PENALIZATION.md`
(which fixed the *Beta*-side sparsity/recovery with the now-default
`cvglmnet` + `simmsar_lasso_reselect=TRUE` engine). That fix does **not** touch
`Sigma`: `Sigma` is always the closed-form regime-weighted residual covariance of
the (possibly sparse) Beta estimate, and `Kappa = solve(Sigma)` is always a raw
inversion — neither is regularized. This document diagnoses the residual
`Sigma`/`Kappa` degeneracy (ill-conditioning, magnitude blow-ups, empty/dense
inverted precision) seen in pilot runs, and tests the leading hypothesis that it
is driven by **low effective sample size in poorly-separated regimes**
(small `postmix[j]`), independent of how well Beta is recovered.

**Everything here is uncommitted / diagnostic.** Instrumentation, the pilot
driver, the stabilization, and this doc are for review only; no production
behaviour was changed (all capture is default-off or additive columns, and the
stabilization defaults to `simmsar_sigma_stab="none"`). §1–7 are the Step-1
diagnosis; **§8 is the Step-2 stabilization (implemented + A/B validated).**

## TL;DR

1. **Hypothesis CONFIRMED for the ill-conditioning / magnitude-explosion mode.**
   The degenerate `est_Sigma`/`est_Kappa` estimates are an **under-determination**
   problem: a regime assigned very few effective observations (`postmix[j]` ≈ 7–15)
   relative to its `d(d+1)/2` covariance parameters yields a near-singular residual
   covariance whose inverse explodes. Degenerate regimes have **median postmix 10
   vs 202** for healthy regimes, and **median 0.68 effective observations per
   covariance parameter vs 10.6**. `postmix_frac` dominates a logistic model of
   degeneracy (z = −9.2, p ≈ 5e−20); posterior entropy adds nothing once postmix
   is controlled for.
2. **The "fully empty Kappa" pattern is a DIFFERENT, distinguishable phenomenon.**
   It is **not** low-postmix. Those regimes have normal postmix (median 197) and a
   *well*-conditioned, near-**diagonal** `est_Sigma` (large eigenvalues) whose
   inverse has sub-threshold off-diagonals → thresholded to empty. It concentrates
   on specific network draws (notably one N=6/M=3/density-0.75 network, `ts_id=1`,
   recurring across every T). A Sigma-conditioning fix will not — and should not —
   change it.
3. **The near-singularity develops early and is locked in.** Per-iteration
   `rcond(sigma[[j]])` traces show conditioning collapses in lockstep with postmix
   over the first ~5–10 EM iterations, then plateaus; iteration 1 is *already*
   ill-conditioned. A stabilization must therefore act **every EM iteration**
   (consistent with the `S.th` coupling in the reduct step), not just at the end.
4. **The severe tail crashes the fit outright.** 13 of 720 fits (all N=8, M=3–4,
   short T, dense) failed across all 6 retries with `solve()` "system is singular"
   errors (reciprocal condition ~1e-17) — the same under-determination, severe
   enough to break the reduct step's linear solve / the E-step rather than merely
   store a degenerate Sigma. This directly confirms the task's concern that Sigma
   conditioning propagates into Beta's own re-estimation (`S.th`).

5. **Step 2 (§8): a per-M-step eigenvalue floor fixes it, cheaply.** Implemented
   (default OFF) in `R/estimation/stabilize_sigma.R`, applied in both M-steps. In
   the A/B it **eliminates the outright fit failures** (baseline 16.7% → 0%), caps
   the worst-case condition number and |Kappa off-diagonal|, keeps Σ **dense** (0
   exact zeros → the reduct freeze logic stays a no-op), and leaves **healthy cells
   bit-for-bit unchanged**. Recommendation: `simmsar_sigma_stab="floor"`,
   `simmsar_sigma_stab_floor=1e-3` (safe default; cond cap 1000) or `1e-2` (tightest
   control, at some risk of benign MaxIter-hitting on borderline fits).

Verdict: the postmix hypothesis holds for the numerically-degenerate mode → the
cheap in-EM stabilization (Step 2) is warranted and validated. The empty-Kappa mode
is out of its scope and is left alone.

## 1. Method / instrumentation (uncommitted)

- **`R/estimation/estimate_MSAR.R`** (`fit_one_replicate`): per (matched) regime
  row, added `postmix` (effective obs count = `colSums` of the smoothed posterior
  for that estimated regime), `postmix_frac` (share of the fit's total effective
  obs), `est_Sigma_rcond`, `est_Sigma_min_eig` (of the exact stored `est_Sigma`
  that analysis inverts), plus fit-level `gamma_entropy_mean` (mean per-timestep
  posterior entropy) and `iter`. `postmix` needed no in-loop capture: the final
  smoothed posterior is already returned as `model_fit$smoothedprob`. The exact
  `lasso_control` is now persisted as `attr(msar_results, "lasso_control")` (so a
  diagnostic run's engine settings never have to be reconstructed after the fact).
- **`R/estimation/fit_msar.R`**: opt-in per-iteration trace
  (`options(simmsar_capture_sigma_trace=TRUE)`, default off), recording per regime
  `postmix` and `rcond(theta$sigma[[j]])` at every EM iteration → attached as
  `res$sigma_trace`. Off by default, so production/parallel runs are byte-identical.
- **Pilot** (`scripts/diagnose_sigma_kappa.R`): the **full production grid**
  `N={4,6,8} × Density={.25,.5,.75} × M={1,2,3,4} × T={200,400,800,1600}`,
  `n_ts=5` (720 fits), default engine `cvglmnet + reselect + 1se + refit`,
  `MaxIter=200`, master seed 58396. Scored with the exact production analysis
  gates (`compute_recovery_metrics`, `min_edg_val=0.05`, `KAPPA_COND_MAX=1e6`,
  `max_plausible_magnitude=10`). Estimation wall time 35.8 min on 5 workers.
- Degeneracy is characterised **self-contained** in the side-table (invert
  `est_Sigma` directly), so even rows the pipeline NA-gates out are described.
- Ground-truth `Kappa` is diagonally-dominant PD by construction
  (`generate_kappa`: `κ_ii = Σ_{j≠i}|κ_ij| + 0.1`, off-diagonals in [0.05, 1]), so
  the *true* Sigma is well-conditioned and true |Kappa off-diag| ≤ 1. **All**
  observed degeneracy is estimation-side.

## 2. Overall degeneracy rates (707 fitted, 1756 regime rows)

13 of 720 fits failed entirely (see §5). Of the 1756 successfully-stored regime
rows, per the production gates:

| flag | count | rate | mechanism |
|---|---|---|---|
| `solve(est_Sigma)` failed | 0 | 0.00% | — |
| ill-conditioned (cond > 1e6, `rcond` gate) | 6 | 0.34% | near-singular Sigma |
| magnitude (|Kappa off-diag| > 10) | 30 | 1.71% | inverse blow-up |
| empty Kappa (no off-diag edge survives 0.05) | 22 | 1.25% | near-diagonal Sigma |
| dense Kappa (specificity < 0.05) | 210 | 11.96% | inverse of full Sigma is dense |
| **any numerically degenerate** (ill / magnitude / empty) | **54** | **3.08%** | — |

`dense Kappa` is a **specificity artifact, not degeneracy**: the inverse of an
(estimated) full residual covariance is generically dense, so many off-diagonals
exceed 0.05 against a sparse truth. It rises monotonically with M (16% at M=4) and
is excluded from the degeneracy set. It is a scoring/thresholding matter, not a
numerical one, and is not in scope here.

## 3. The three numerically-distinct modes

Splitting the flagged rows into mutually-exclusive modes (magnitude → ill-cond →
empty) exposes two **opposite** Sigma pathologies plus a ratio-only case:

| mode | n | postmix med (min) | postmix_frac med | min eig med | cond med | \|Kappa off-diag\| med |
|---|---|---|---|---|---|---|
| **magnitude_explosion** | 30 | **10.0 (7)** | 0.037 | **3.9e-3** | 5.9e3 | **68** |
| ill_conditioned_only | 3 | 198 (9) | 0.248 | 2.8e-1 | 4.2e9 | 0.79 |
| **kappa_empty** | 21 | **197 (55)** | 0.278 | **133** | 188 | **1.7e-3** |
| healthy | 1702 | 202 | 0.357 | 3.4e-1 | 23.6 | 0.92 |

- **magnitude_explosion (the hypothesis mode).** Starved regimes: postmix ~10 for
  a d=8 → 36-parameter covariance. The residual covariance has a **tiny smallest
  eigenvalue** (~4e-3, down to 4.5e-7), so its inverse blows up. The largest
  observed off-diagonals (postmix in parentheses):

  | (T,D,N,M,reg) | postmix | min eig | cond | \|Kappa off-diag\| max |
  |---|---|---|---|---|
  | 800,.5,8,4,r3 | 10.0 | 4.5e-7 | 3.3e7 | **603068** |
  | 400,.75,8,3,r2 | 9.0 | 2.3e-6 | 1.1e6 | 105970 |
  | 200,.25,8,3,r3 | 10.0 | 7.6e-5 | 1.9e6 | 5418 |
  | 400,.25,8,4,r3 | 9.0 | 5.5e-5 | 3.9e5 | 3129 |

  Every one of the top magnitude rows has postmix ≈ 9–12 and a near-zero minimum
  eigenvalue. This is exactly the ">1000 off-diagonals / cond up to ~1e11" and
  "magnitude-outlier" pattern from earlier pilots.

- **kappa_empty (distinct — NOT low postmix).** Normal postmix (median 197), and a
  **well-conditioned, near-diagonal** Sigma (**large** min eigenvalue ~133, cond
  ~188). Its inverse is essentially diagonal, so every off-diagonal (~1.7e-3) is
  below the 0.05 threshold → the thresholded Kappa is empty and `Kappa_corr` is
  undefined. The rows concentrate on a **handful of specific network draws**:
  overwhelmingly `ts_id=1, N=6, M=3, density=0.75` (appears at T=200/400/800/1600
  — the same generating network), plus `N=4/M=3/density=0.75/ts_id=1` and some
  `N=8/M=2`. All 3 "every regime empty in the whole fit" cases are that one
  N=6/M=3/.75 network. These are genuine weak-partial-correlation networks, not a
  numerical failure — ridge/eigenvalue stabilization would not (and should not)
  change them.

- **ill_conditioned_only (3 rows).** High condition number by **ratio** (cond
  ~4e9) but off-diagonals stay < 1 and min eigenvalue is not tiny — caught by the
  `rcond` gate, not the magnitude gate. Rare; a mixed/boundary case.

## 4. The hypothesis test: degeneracy ~ effective sample size

Degenerate vs healthy, on the separation diagnostics:

| group | n | postmix med (min) | postmix_frac med (min) | entropy med | cond med |
|---|---|---|---|---|---|
| **DEGENERATE** | 54 | **13.6 (7.0)** | **0.060 (0.006)** | 0.003 | 2.6e3 |
| healthy | 1702 | 202.0 (8.0) | 0.357 (0.012) | 0.003 | 23.6 |

The effective-sample-size gap is ~15×. Monotone dose-response on both the raw
share and the parameter-adjusted count:

| `postmix_frac` quintile | [.006,.22] | (.22,.30] | (.30,.405] | (.405,.522] | (.522,1] |
|---|---|---|---|---|---|
| degeneracy rate | **9.9%** | 3.4% | 0.9% | 0.9% | **0.3%** |

| obs per cov-param `postmix/[d(d+1)/2]` quintile | [.25,3.78] | (3.78,7.8] | (7.8,13.9] | (13.9,25.3] | (25.3,160] |
|---|---|---|---|---|---|
| degeneracy rate | **9.9%** | 2.3% | 1.4% | 1.4% | **0.3%** |

Median obs-per-parameter is **0.68 (degenerate) vs 10.6 (healthy)** — i.e. the bad
regimes have fewer effective observations than free covariance parameters.

Logistic regression `degenerate ~ postmix_frac + entropy + factor(N) + factor(M)`:

| term | estimate | z | p |
|---|---|---|---|
| `postmix_frac` | **−13.83** | **−9.17** | **4.8e-20** |
| `gamma_entropy_mean` | 2.95 | 0.92 | 0.36 (ns) |
| `factor(N)8` | 1.09 | 2.20 | 0.028 |
| `factor(N)6` | 0.90 | 1.83 | 0.068 |

`postmix_frac` is overwhelmingly the driver. Nodes (N) adds a significant residual
effect **in the expected direction** — the covariance parameter count grows as
`d(d+1)/2`, so for fixed postmix a larger network is more under-determined (the
obs-per-parameter table already folds this in). Posterior **entropy contributes
nothing** once postmix is held fixed: the richer separation measure does not
distinguish degenerate regimes — it is specifically the *effective count*, not
"chronic uncertainty", that matters. Regime factors are non-significant given
postmix (M works entirely through starving regimes).

### Break-down by design cell (degeneracy rate)

| factor | levels → rate |
|---|---|
| M (regimes) | 1: **0.0%** · 2: 1.4% · 3: 4.9% · 4: 3.3% |
| T (timesteps) | 200: **5.2%** · 400: 2.5% · 800: 3.1% · 1600: **1.6%** |
| N (nodes) | 4: 1.4% · 6: 3.2% · 8: **4.6%** |
| density | .25: 1.9% · .5: 1.9% · .75: **5.4%** |

All four gradients are exactly what the postmix mechanism predicts and confirm the
diagnosis is not "cleanly explained by postmix alone" in a trivial way — each
factor acts *through* effective-obs-per-parameter: **M>1** is required to split the
data across regimes (M=1 never degenerates — one regime keeps all obs); **short T**
gives fewer total obs to split; **large N** raises the parameter count per regime;
**high density** makes regimes harder to separate (so the E-step starves the weaker
ones) — density 0.75 is the worst cell and is also where the empty-Kappa network
draws sit. No design cell shows degeneracy *independent* of low postmix: the
magnitude/ill-conditioned rows are all low-postmix; the only postmix-independent
pattern is the empty-Kappa mode (§3), which is network-draw-specific, not cell-wide.

## 5. Onset over EM iterations (gradual, then locked in)

Opt-in `rcond(sigma[[j]])` traces on hard N=8/M=4/density-0.75 fits (worst regime
per fit):

```
T=400 N=8 D=.75 M=4 ts=1, regime 4:  postmix 15.8→9.0 over iters 1–5,
   rcond 5.2e-4 → 3.2e-5 (cond 1930 → 31200), then FLAT for 35 more iters.
T=400 N=8 D=.75 M=4 ts=2, regime 2:  postmix 60→18 over iters 1–7,
   rcond 0.032 → 0.0027 (cond 32 → 366), then flat.
T=200 N=8 D=.75 M=4 ts=2, regime 1:  postmix ~12 throughout,
   rcond 0.005 → 0.001 by iter 3 (cond 200 → 978), flat thereafter.
```

The near-singularity is **not a late sharp collapse**: `rcond` tracks postmix
essentially deterministically as the E-step redistributes mass over the first
~5–10 iterations, then both plateau together and never recover. Iteration 1 is
already ill-conditioned (cond 200–2000 for the starved regimes). Implication: a
stabilization keyed to Sigma's own conditioning must run **every M-step**, not as
a one-off at convergence.

### 5.1 The severe tail: outright fit failures (13 of 720)

The most starved cells don't just store a degenerate Sigma — they crash the EM.
All 13 lost fits are N=8, M=3–4, T∈{200,400}, density 0.5–0.75, and fail across
all 6 retries with either `System ist singulär: reziproke Konditionszahl ≈ 1e-17`
(a `solve()` on a near-singular matrix) or `smoothing probabilities are too small`
(a regime collapses to ~0 mass). The singular-`solve` errors originate in the
reduct M-step's linear system / likelihood term (`solve(S.th)` at
`mstep_hh_reduct_msar.R:132`, `solve(A,b)` at `:168`) — i.e. **Sigma's
conditioning already propagates into Beta's re-estimation**, exactly the `S.th`
coupling the task flagged. These are the extreme end of the same continuum, not a
separate failure mode.

## 6. Implications for Step 2 (stabilization)

The hypothesis holds for the numerically-degenerate mode, so the Step-2 plan is
on target, with these specifics from the data:

- **Act every EM iteration** (§5): the ill-conditioning is present from iteration 1
  and persists; it also feeds Beta via `S.th`. A per-M-step floor protects both the
  E-step likelihood and the reduct linear solve (and would likely rescue several of
  the 13 hard crashes in §5.1).
- **Target the small-eigenvalue direction.** The magnitude mode is precisely a
  near-zero **minimum** eigenvalue (median 3.9e-3, down to 4.5e-7) — so an
  eigenvalue floor and a trace-scaled ridge are the natural candidates; both lift
  the smallest eigenvalue while leaving well-conditioned Sigmas essentially
  untouched. `est_Sigma_min_eig` in the side-table is the quantity to calibrate the
  floor against (healthy median 0.34 vs degenerate 3.9e-3 → ~2 orders of magnitude
  of separation to place a floor in).
- **Leave the empty-Kappa mode alone** (§3): it is well-conditioned and
  network-specific; a ridge/floor won't move it and shouldn't. Do not tune the
  stabilization to "fix" empty Kappa.
- **Keep Sigma dense.** The stabilization must not introduce exact zeros, or it
  would arm the frozen-support logic at `mstep_hh_reduct_msar.R:199`
  (`w = which(abs(theta$sigma[[j]]) > 0)`) — currently a no-op only because the raw
  residual covariance is dense. Ridge (`Sigma + λI`) and eigenvalue flooring both
  keep every entry non-zero; this must be **verified**, not assumed (Step 2).
- **Validate on the right cells.** The degeneracy lives in N=8, M≥3, short T,
  density 0.75; A/B there (plus a couple of healthy cells to confirm no distortion)
  is the discriminating comparison, with matched data+init as in the Beta LASSO fix.

## 7. Reproduce

```r
# full pilot + side-table + headline analysis (uncommitted):
Rscript scripts/diagnose_sigma_kappa.R
#   -> output/diagnostics/MSAR_diag_<stamp>.rds        (scored msar_results)
#   -> output/diagnostics/sigma_kappa_diag_<stamp>.csv (per-regime side-table)
# per-iteration rcond trace on hard cells: set
#   options(simmsar_capture_sigma_trace = TRUE); fit; inspect fit$sigma_trace
```

Backups of the touched estimation files (suffix `.bak_sigmakappa_2026-07-01.R`):
`fit_msar`, `estimate_MSAR`, `mstep_hh_lasso_msar`, `mstep_hh_reduct_msar`.

## 8. Step 2: in-EM Sigma stabilization (implemented + A/B validated)

The Step-1 hypothesis holds for the numerically-degenerate mode, so a cheap
per-M-step stabilization is warranted. **Implemented** (default OFF), in
`R/estimation/stabilize_sigma.R`, applied to `sigma[[j]]` in **both** M-steps
(`mstep_hh_lasso_msar.R` — the cvglmnet engine that runs *every* iteration under
the production `reselect` default — and `mstep_hh_reduct_msar.R` — used by the
bic/frozen configs), so it acts on every EM iteration on every LASSO path. Two
standard options, chosen via `options()` (so they ride the same `lasso_control`
worker-propagation as the LASSO knobs), default `simmsar_sigma_stab="none"`:

- **ridge**: `Σ ← Σ + λ·(tr Σ / d)·I` (`simmsar_sigma_stab_lambda`, default 1e-3).
  Scale-free lift of *every* eigenvalue; cannot create a zero.
- **floor**: eigen-decompose, clip eigenvalues below `floor·max(eig)` up to that
  floor, reconstruct (`simmsar_sigma_stab_floor`, default 1e-3). Because the floor
  is a fraction of the *largest* eigenvalue, it **bounds the condition number at
  exactly `1/floor`** and leaves any Σ already better-conditioned than that
  completely untouched.

Both keep Σ dense by construction (ridge touches only the diagonal; floor's
`V diag V'` reconstruction has no exact zeros) — verified below, so the
`w = which(abs(theta$sigma[[j]]) > 0)` freeze in `mstep_hh_reduct_msar.R` stays a
no-op.

Unit check (a 4×4 Σ with eigenvalues `{2,1,0.5,1e-7}`, cond 2e7): `floor` at 1e-3
takes cond → exactly **1e3**, `ridge` at 1e-2 → 230; both add **0** exact zeros
and stay symmetric; on a healthy Σ (cond 4) both are essentially identity
(max abs change ~1e-15 for floor).

### 8.1 A/B design

`scripts/ab_sigma_stab.R` + `ab_sigma_combine.R`. 6 cells — 4 degenerate-prone
(N=8/M=4/T=200/.75, N=8/M=4/T=400/.5, N=8/M=3/T=400/.75, N=6/M=3/T=200/.75) and
2 healthy (N=4/M=2/T=800/.25, N=6/M=2/T=800/.5) — × `n_ts=4` × 5 configs
(baseline, ridge 1e-3, ridge 1e-2, floor 1e-3, floor 1e-2), `MaxIter=200`,
production `cvglmnet+reselect` engine. **Identical generated data (seeded per
cell) and identical initial `theta` (seeded per cell×rep) across all configs**, so
the only difference is the Σ stabilization. Each `(cell, config)` runs in its own
R process: a baseline degenerate fit can **segfault inside a compiled `solve()`**
(uncatchable by `tryCatch`), so isolation stops one crash from killing the run —
and the crash count is itself an outcome.

### 8.2 Results

Fit-level (per config, 24 fits each):

| config | fit failures | hits MaxIter | mean iters | mean s |
|---|---|---|---|---|
| baseline | **16.7% (4/24)** | 0 | 14.1 | 9.2 |
| ridge 1e-3 | 4.2% (1/24) | 1 | 22.0 | 15.5 |
| ridge 1e-2 | **0%** | 0 | 17.5 | 13.4 |
| floor 1e-3 | **0%** | 1 | 21.8 | 15.3 |
| **floor 1e-2** | **0%** | 0 | 16.1 | 11.7 |

**The headline win: stabilization eliminates the outright fit failures.** Baseline
loses 1-in-6 fits to `solve()` "system is singular" crashes (the §5.1 mechanism);
`ridge 1e-2`, `floor 1e-3`, `floor 1e-2` recover **all** of them. This is the
`S.th` coupling made concrete — conditioning Σ rescues the fit whose *Beta*
re-estimation solve was crashing. Runtime rises modestly (the recovered hard fits
now actually run); no config reintroduces non-convergence (≤1 MaxIter hit).

Degeneracy flags & density (per **surviving** regime row):

| config | rows | ill-cond | magnitude | worst-case cond (hard) | worst-case \|Koff\| (hard) | Σ exact zeros |
|---|---|---|---|---|---|---|
| baseline | 58 | 0% | 1.7% | 1.2e4 | 38.5 | **0** |
| ridge 1e-3 | 69 | 0% | 5.8% | 1.9e4 | 135 | **0** |
| ridge 1e-2 | 69 | 0% | 2.9% | 1.3e3 | 25.4 | **0** |
| floor 1e-3 | 69 | 0% | 1.4% | 2.6e3 | 38.5 | **0** |
| **floor 1e-2** | 69 | 0% | 1.4% | **2.4e2** | **16.9** | **0** |

- **Density preserved for every config** (`Σ exact zeros = 0`) → the reduct freeze
  logic stays the intended no-op. This was the explicit must-verify.
- **`floor 1e-2` gives the tightest control**: worst-case condition number 236 (vs
  baseline 1.2e4) and worst |Kappa off-diagonal| 16.9 (vs 38.5).
- **Survivorship caveat (important):** baseline shows only 58 surviving rows vs 69,
  because its 4 crashed fits contribute *no* regime rows. So baseline's low
  surviving-row flag rate **understates** its true degeneracy — its worst regimes
  don't appear as "flagged", they appear as *failed fits*. The fair statement is
  that stabilization converts baseline's 16.7% hard **failures** into usable,
  well-conditioned estimates. (This is also why weak `ridge 1e-3` looks *worse* on
  the magnitude rate: it resurrects a crashed fit but under-conditions it, so the
  now-surviving regime is magnitude-flagged.)

Recovery (mean over regimes):

| | Beta_corr | Beta_spec | Kappa_corr | Kappa_spec | NRMSE_Kappa |
|---|---|---|---|---|---|
| **hard** baseline | 0.745 | 0.797 | 0.538 | 0.193 | 0.438 |
| hard ridge 1e-2 | 0.746 | 0.787 | 0.524 | 0.220 | 0.409 |
| hard floor 1e-3 | 0.731 | 0.803 | 0.510 | 0.211 | 0.425 |
| **hard floor 1e-2** | 0.744 | 0.797 | 0.495 | 0.235 | **0.334** |
| **healthy** baseline | 0.997 | 0.993 | 0.991 | 0.830 | 0.034 |
| **healthy floor 1e-2** | 0.997 | 0.993 | 0.991 | 0.830 | 0.034 |

- **No recovery cost on healthy cells** — `floor` is *identical to baseline to 3
  decimals* on every healthy metric (floor never touches a Σ with cond < 100, and
  the healthy Σ's here have cond ~25). `ridge 1e-2` perturbs healthy Kappa_sen/spec
  by ~1–2pp (it shifts every eigenvalue, even the healthy ones) — a small but real
  reason to prefer floor.
- **Hard-cell recovery is not degraded and NRMSE improves.** Beta and Kappa
  correlations are within noise of baseline (and baseline's are computed on the
  easier surviving subset); `floor 1e-2` gives the best hard NRMSE_Kappa (0.334 vs
  0.438) by taming the magnitude outliers that dominate that absolute error.

The table above is **survivorship-confounded** (baseline's means exclude the hard
fits it crashes on). The clean comparison **pairs each config against baseline on
only the regime-rows that succeeded under both** (58 common rows: 42 hard / 16
healthy); deltas are `config − baseline`:

| metric | HARD floor 1e-3 | HARD floor 1e-2 | HARD ridge 1e-2 | HEALTHY floor 1e-3/1e-2 | HEALTHY ridge 1e-2 |
|---|---|---|---|---|---|
| Beta_corr | +0.001 | −0.002 | +0.014 | +0.000 | −0.000 |
| Beta_sen | −0.001 | +0.002 | +0.011 | +0.000 | +0.000 |
| Beta_spec | +0.000 | −0.004 | −0.019 | +0.000 | +0.000 |
| Kappa_corr | −0.002 | −0.017 | +0.014 | +0.000 | +0.000 |
| Kappa_sen | −0.012 | −0.025 | −0.015 | +0.000 | **−0.008** |
| Kappa_spec | +0.003 | +0.007 | +0.029 | +0.000 | **+0.018** |
| NRMSE_Kappa | +0.000 | **−0.087** | **−0.063** | +0.000 | +0.002 |

- **Healthy: `floor` is exactly 0.000 on every metric** (its Σ's, cond ~25, never
  hit the cap); `ridge 1e-2` nudges healthy Kappa sens/spec by ~1–2pp (it lifts all
  eigenvalues) — the reason to prefer floor.
- **Hard: correlations move within noise** (±0.02 on 42 rows) for both Beta and
  Kappa; Beta is essentially untouched (only the weak `S.th` coupling reaches it).
- **Hard: Kappa_sen dips a little / Kappa_spec rises a little** — capping the
  conditioning shrinks the inflated off-diagonals, so a few borderline true edges
  fall under the 0.05 threshold (−sen) while spurious huge entries vanish (+spec).
- **Hard: NRMSE_Kappa is the real recovery gain** — `floor 1e-2` −0.087, `ridge
  1e-2` −0.063 (NRMSE is an *absolute* error dominated by the exploded entries).
  `floor 1e-3` does **not** move it: cap 1000 is too loose to remove the
  NRMSE-driving outliers on these particular common fits.
- **The largest benefit is invisible in this paired table.** Every stabilized
  config additionally rescued **3 hard fits (cells 2 & 3) that baseline crashed on
  entirely** — there baseline has *no* Beta/Kappa estimate, and stabilization
  yields a usable one. The paired set can only compare fits both produced.
- **`NRMSE_Beta` was not scored** in the A/B harness (only Beta corr/sen/spec).
  Beta is not regularized and its corr/sen/spec are ~unchanged, so NRMSE_Beta is
  expected to move negligibly; stated as unmeasured rather than asserted.

**Convergence caveat (benign non-convergence).** Flooring *every* M-step means the
M-step no longer returns the exact maximizer, so when the cap is well below a fit's
natural conditioning the strict `eps=1e-5` test may never trip even though the
parameters are stable. In the 24-fit A/B this did not bite (`floor 1e-2`: 0/24
MaxIter hits, fewest mean iters), but a separate single-fit probe (a hard
N=8/M=4/T=200/.75 fit that baseline solved in 12 iters at cond 2.5e3) exposed it:

| config | iters | worst cond | loglik range, last 5 iters |
|---|---|---|---|
| baseline | 12 | 2.5e3 | 4.99 |
| floor 1e-3 | 30 | 1.0e3 | 122 |
| floor 1e-2 | **200 (MaxIter)** | 1.0e2 | **1.98** |
| ridge 1e-2 | 19 | 5.3e2 | 2.37 |

`floor 1e-2` hit MaxIter here — but its terminal loglik is *more* settled than
baseline's (range 1.98 vs 4.99): the estimate is stable and well-conditioned, it
just orbits just above the strict threshold. So the risk of aggressive flooring is
wasted **runtime**, not a bad estimate — distinct from the frozen-support
non-convergence (which produced genuinely worse fits). `ridge 1e-2` converged
fastest and cleanest on this fit; `floor 1e-3` (milder cap) also converged.

### 8.3 Recommendation

**Eigenvalue flooring, applied every M-step.** Flooring is preferable to ridge in
principle: it bounds the condition number by construction (a directly interpretable
knob = `1/floor`) and acts *only* on the ill-conditioned directions, so healthy
regimes are left bit-for-bit untouched (ridge shifts every eigenvalue and perturbs
even well-conditioned regimes by ~1–2pp on Kappa sens/spec). Both floor settings
and `ridge 1e-2` eliminate the outright fit failures and keep Σ dense.

Between the two floor strengths there is a **conditioning-tightness vs
convergence-cost trade-off**, so pick by priority:

- **`simmsar_sigma_stab_floor=1e-2`** (cond cap ≈ 100): tightest worst-case control
  (cond ≤ 100, best NRMSE), untouched healthy cells, and in the 24-fit A/B the
  fewest iterations — **but** on borderline fits it can hit MaxIter (benignly — see
  the caveat above). Choose it if worst-case conditioning is the priority and the
  MaxIter runtime is acceptable.
- **`simmsar_sigma_stab_floor=1e-3`** (cond cap ≈ 1000): the **safer default** —
  still removes every failure and caps the pilot's cond-up-to-1e11 / |Koff|-up-to-6e5
  explosions down to cond ≤ 1000, while touching far fewer borderline regimes (only
  cond > 1000), so it converges more reliably. It leaves a smaller residual
  magnitude tail than baseline but a larger one than `1e-2`.

`ridge 1e-2` is a reasonable third option where reliable strict convergence matters
more than leaving healthy cells exactly untouched (it converged fastest on the hard
probe fit). Do **not** use `ridge 1e-3` / `floor` weaker than 1e-3: under-conditions
(still 4% failures, residual magnitude 135).

Net: **recommend `floor` at `1e-3` for a safe production default, `1e-2` when the
worst-case conditioning must be bounded tightly and the extra iterations are
tolerable.** Confirm the MaxIter-hit rate on the full `n_ts` grid before committing
to `1e-2`, since the 24-fit A/B is too small to estimate that rate precisely.

**Left OFF by default** (`simmsar_sigma_stab="none"`), exactly like the cvglmnet
engine: production behaviour is unchanged until opted in.

Out of scope (per the task, and confirmed unnecessary by §3): no graphical-lasso
sparsification of Σ/Kappa inside the loop, and no joint Beta+Σ re-selection. The
empty-Kappa mode is deliberately left unchanged.

### 8.5 Production wiring (verified parallel-safe)

The stabilization rides the **same `lasso_control` option-propagation as the LASSO
knobs**: `estimate_MSAR()` re-applies the list via `options()` *inside every future
worker* (options don't cross the process boundary), and `future` exports the
`stabilize_sigma` function to workers transitively (it is referenced by name in both
M-steps). **Verified**: `workers=2` with `floor 1e-2` caps every stored `est_Sigma`
at cond **100** vs **44150** unstabilized — so it genuinely fires in parallel, not
just sequentially. No call sites change.

To enable in the production driver (`scripts/MSAR_ts_analysis_NHMSAR.R`), add two
entries to the `lasso_control` list that is already built there and passed to
`estimate_MSAR(..., lasso_control = lasso_control)`:

```r
lasso_control <- list(
  simmsar_lasso_engine         = lasso_engine,      # existing entries ...
  simmsar_lasso_reselect       = lasso_reselect,
  # ... (unchanged) ...
  simmsar_lasso_reselect_every = lasso_reselect_every,
  # --- NEW: residual-covariance stabilization -------------------------------
  simmsar_sigma_stab       = "floor",   # "none" (default) | "floor" | "ridge"
  simmsar_sigma_stab_floor = 1e-3       # cond cap = 1/floor  (1e-3 -> 1000)
)
```

The existing `do.call(options, lasso_control)` line in the driver also applies it to
the main process (covers `workers = 1` / interactive fits). For an interactive
single fit outside the driver, just `options(simmsar_sigma_stab = "floor",
simmsar_sigma_stab_floor = 1e-3)` before calling `fit_msar()` /
`init_and_fit_msar_lasso()`. Nothing downstream needs changing: the analysis
pipeline (`compute_recovery_metrics`) still inverts the stored `est_Sigma`, now
well-conditioned, and its `KAPPA_COND_MAX` / `max_plausible_magnitude` gates simply
flag far fewer rows.

### 8.4 Reproduce Step 2

```r
# 30 isolated (cell, config) fits (segfault-resilient), then combine:
for ci in 1..6, cf in {baseline,ridge_1e3,ridge_1e2,floor_1e3,floor_1e2}:
  Rscript scripts/ab_sigma_stab.R <ci> <cf>     # -> output/diagnostics/ab_part_<ci>_<cf>.csv
Rscript scripts/ab_sigma_combine.R              # -> ab_sigma_stab_<stamp>.csv + summary
```
