# Reproducibility Guide for SimMSAR

## Quick Answer

✅ **Yes, setting `set.seed()` once at the start of `MSAR_ts_analysis_NHMSAR.R` is sufficient!**

I've uncommented it for you (line 99).

---

## Where Randomness Occurs

Your pipeline has random operations at multiple stages:

### 1. Network Generation (`generate_Beta.R`, `generate_kappa.R`, etc.)
```r
# Random edge placement
indices <- sample(N*N, num_nonzero, replace = FALSE)

# Random edge weights
values <- runif(num_nonzero, min_edg_val, max_edg_val) * signs
```

### 2. Transition Matrix (`generate_transmat.R`)
```r
# Random transition probabilities
remain_prob <- runif(1, min = remain_lower, max = remain_upper)
```

### 3. Time Series Simulation (`generate_timeseries.R`)
```r
# Initial values
init <- runif(N[j], min = 0, max = 5)

# Regime switching
reg_index <- sample.int(M[k], 1, prob = transmat[reg_index, ])

# Gaussian noise
mvtnorm::rmvnorm(1, rep(0, N[j]), curreg_sigma)
```

### 4. Model Initialization (`init_theta_msar.R`)
```r
# k-means initialization
centers <- d.data[sample(1:dim(d.data)[1], M), ]

# Random prior and transition matrix
prior <- NHMSAR:::normalise(runif(M))
transmat <- NHMSAR:::mk_stochastic(diag(1,M) + matrix(runif(M*M), M))
```

---

## Why One Seed Works

R's random number generator maintains a **global state**. When you call `set.seed()`:

1. ✅ All subsequent `runif()`, `rnorm()`, `sample()`, etc. calls use this state
2. ✅ The state updates deterministically with each random call
3. ✅ Same seed → same sequence of random numbers → same results

### Example:
```r
set.seed(123)
runif(3)  # [1] 0.2875775 0.7883051 0.4089769

set.seed(123)  # Reset to same seed
runif(3)  # [1] 0.2875775 0.7883051 0.4089769  # Same results!
```

---

## When One Seed Is NOT Enough

You need separate seeds or special handling if you have:

### ❌ Parallel Processing
```r
# DON'T do this with simple set.seed()
library(parallel)
mclapply(1:4, function(i) {
  rnorm(10)  # Each worker might have same seed!
})

# DO this instead
library(parallel)
RNGkind("L'Ecuyer-CMRG")  # Use parallel-safe RNG
set.seed(123)
mclapply(1:4, function(i) {
  rnorm(10)  # Now properly independent
}, mc.cores = 4)
```

### ❌ Distributed Computing (future, foreach)
```r
# Need special setup for distributed RNG
library(doRNG)
registerDoRNG(123)  # Instead of set.seed()
```

### ✅ **Your code has NONE of these!**
You use sequential loops everywhere, so one seed is perfect.

---

## Best Practices for Your Workflow

### 1. **Current Setup (Recommended)**
```r
# In MSAR_ts_analysis_NHMSAR.R
set.seed(83742)  # Line 99

# Then run entire pipeline
Timeseries_data <- generate_timeseries(...)
MSAR_dynamics_list <- estimate_MSAR(...)
```

✅ **This gives you:**
- Reproducible network structures
- Reproducible time series
- Reproducible model initialization
- Reproducible EM convergence (if stochastic components exist)

### 2. **For Experiments with Different Conditions**

If you want to test different parameters but keep randomness consistent:

```r
# Same random structures, different parameters
set.seed(83742)
results_short <- estimate_MSAR(T = 1000, ...)

set.seed(83742)  # Reset seed
results_long <- estimate_MSAR(T = 5000, ...)
# Now networks are identical, only T differs
```

### 3. **For Multiple Independent Runs**

If you want statistically independent replications:

```r
seeds <- c(123, 456, 789, 101112)

results_list <- list()
for (i in seq_along(seeds)) {
  set.seed(seeds[i])
  results_list[[i]] <- generate_timeseries(...)
}
# Each run is reproducible but independent
```

---

## Verifying Reproducibility

Test that your seed works:

```r
# Run 1
set.seed(83742)
ts1 <- generate_timeseries(Density = 0.3, N = 4, M = 2, T = 100, n_ts = 1, ...)
first_value_run1 <- ts1$timeseries_data[[1]][1,1]

# Run 2 (same seed)
set.seed(83742)
ts2 <- generate_timeseries(Density = 0.3, N = 4, M = 2, T = 100, n_ts = 1, ...)
first_value_run2 <- ts2$timeseries_data[[1]][1,1]

# Should be TRUE
identical(first_value_run1, first_value_run2)

# Should be TRUE (entire tibble identical)
all.equal(ts1, ts2)
```

---

## What the Seed Controls

| Component | Reproducible? | Why |
|-----------|--------------|-----|
| Network structure (edges) | ✅ Yes | `sample()`, `runif()` |
| Edge weights | ✅ Yes | `runif()` |
| Transition matrices | ✅ Yes | `runif()` |
| Initial time series values | ✅ Yes | `runif()` |
| Regime sequences | ✅ Yes | `sample.int()` |
| Noise terms | ✅ Yes | `rmvnorm()` |
| k-means initialization | ✅ Yes | `sample()` |
| EM initialization | ✅ Yes | `runif()` |
| EM convergence path | ⚠️ Mostly | Usually deterministic after init |
| Final parameter estimates | ✅ Yes | If EM is deterministic |

---

## Common Issues

### Issue 1: "I set the seed but got different results"

**Possible causes:**
1. ❌ Different R version (RNG changed between versions)
2. ❌ Different package versions (especially mvtnorm, MASS)
3. ❌ Seed set AFTER random operations started
4. ❌ Parallel processing without proper RNG setup

**Solution:**
```r
# At the very top of your script
set.seed(83742)
RNGkind(sample.kind = "Rounding")  # For R >= 3.6.0 backward compatibility
```

### Issue 2: "Different runs of stat_analysis.R give different results"

**Check:**
- Are there random operations in stat_analysis.R? (permutations, bootstrap?)
- If yes, add `set.seed()` at the start of that script too

---

## Summary

✅ **For your current setup:**
- One `set.seed(83742)` at line 99 of `MSAR_ts_analysis_NHMSAR.R` is **sufficient**
- This makes **everything** reproducible
- No parallel processing → no special RNG needed
- You're good to go!

✅ **To guarantee reproducibility:**
- Use the same R version
- Use the same package versions
- Set seed before running pipeline
- Don't manually modify RNG state

✅ **Your results will be reproducible if:**
- Same seed
- Same parameters
- Same R/package versions
- Sequential execution (no parallel changes)

---

## Session Info Best Practice

Save your session info with results:

```r
# At the end of MSAR_ts_analysis_NHMSAR.R
if (save_output) {
  session_file <- file.path(output_dir, sprintf("sessionInfo_%s.txt", timestamp))
  writeLines(capture.output(sessionInfo()), session_file)
  cat(sprintf("✓ Saved session info: %s\n", session_file))
}
```

This lets you verify:
- R version
- Package versions
- Platform details
- RNG kind

All critical for reproducibility!
