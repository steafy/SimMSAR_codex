# Performance Analysis - SimMSAR

## 🔍 Identified Bottlenecks

### **Critical Issue #1: Growing Vectors in Loops (O(n²) complexity)**

**Location**: `R/estimation/estimate_MSAR.R` lines 328-465

**Problem**: Multiple vectors are grown using `c()` inside nested loops:

```r
# BEFORE (SLOW - O(n²)):
vec_all_Wtemp_cor <- numeric()
vec_all_Wtemp_sen <- numeric()
vec_all_Wcont_cor <- numeric()
# ... 10 total vectors

for (l in 1:n_ts) {
  for (m in 1:M[k]) {
    vec_all_Wtemp_cor <- c(vec_all_Wtemp_cor, new_value)  # COPIES ENTIRE VECTOR!
    vec_all_Wtemp_sen <- c(vec_all_Wtemp_sen, new_value)
    vec_all_Wcont_cor <- c(vec_all_Wcont_cor, new_value)
    # ... etc
  }
}
```

**Why it's slow**: Every `c()` call creates a NEW vector and copies all existing elements. For 30 time series × 2 regimes = 60 iterations, this means:
- 1st iteration: copy 0 elements
- 2nd iteration: copy 1 element
- 3rd iteration: copy 2 elements
- ...
- 60th iteration: copy 59 elements
- **Total copies**: 0+1+2+...+59 = 1,770 unnecessary copies!

**Impact**: With n_ts=30, M=2, this creates ~1,770 vector copies. With n_ts=100, M=4, it's ~79,800 copies!

---

### **Critical Issue #2: Growing Regime Sequence in Loop**

**Location**: `R/generation/generate_timeseries.R` line 199

**Problem**:
```r
# BEFORE (SLOW):
Rseq <- numeric()
for (m in 2:totTime[t]) {
  reg_index <- sample.int(M[k], 1, prob = transmat[reg_index, ])
  Rseq <- c(Rseq, reg_index)  # COPIES ENTIRE VECTOR!
}
```

**Impact**: For T=3500 timesteps, this creates 3,500 vector copies (6+ million element copies total!).

---

### **Moderate Issue #3: Growing Lists in Loops**

**Locations**:
- `R/analysis/stat_analysis.R` line 25: `append(corr_list, list(temp))`
- `R/visualization/Create_3D_surfaceplots_for_means.R` line 30: `append(data_list, list(temp))`

**Problem**: Same issue as above, but with lists instead of vectors.

---

## 🚀 Optimizations

### **Fix #1: Pre-allocate Vectors in estimate_MSAR.R**

```r
# AFTER (FAST - O(n)):
# Pre-allocate to maximum possible size
max_size <- n_ts * M[k]
vec_all_Wtemp_cor <- numeric(max_size)
vec_all_Wtemp_sen <- numeric(max_size)
vec_all_Wcont_cor <- numeric(max_size)
# ... etc

idx <- 0  # Index counter
for (l in 1:n_ts) {
  for (m in 1:M[k]) {
    idx <- idx + 1
    vec_all_Wtemp_cor[idx] <- new_value  # DIRECT ASSIGNMENT!
    vec_all_Wtemp_sen[idx] <- new_value
    vec_all_Wcont_cor[idx] <- new_value
  }
}

# Trim to actual size (removes NAs from failed fits)
vec_all_Wtemp_cor <- vec_all_Wtemp_cor[1:idx]
vec_all_Wtemp_sen <- vec_all_Wtemp_sen[1:idx]
# ... etc
```

**Expected speedup**: **10-50x faster** for this section

---

### **Fix #2: Pre-allocate Regime Sequence**

```r
# AFTER (FAST):
Rseq <- integer(totTime[t])  # Pre-allocate
for (m in 2:totTime[t]) {
  reg_index <- sample.int(M[k], 1, prob = transmat[reg_index, ])
  Rseq[m] <- reg_index  # DIRECT ASSIGNMENT!
}
Rseq <- Rseq[-1]  # Remove first element (initialized to 0)
```

**Expected speedup**: **50-100x faster** for this section

---

### **Fix #3: Pre-allocate Lists**

```r
# AFTER:
corr_list <- vector("list", length = estimated_size)
idx <- 0
for (...) {
  idx <- idx + 1
  corr_list[[idx]] <- temp
}
corr_list <- corr_list[1:idx]  # Trim
```

---

## 📊 Expected Overall Impact

| Component | Current | Optimized | Speedup |
|-----------|---------|-----------|---------|
| estimate_MSAR (vector growth) | ~slow | ~fast | **10-50x** |
| generate_timeseries (Rseq) | ~slow | ~fast | **50-100x** |
| Overall pipeline | baseline | improved | **5-20x** |

**Note**: Actual speedup depends on n_ts and T values. Larger values = bigger speedup.

---

## 🎯 Implementation Priority

1. **HIGH**: Fix estimate_MSAR.R vector growth (biggest impact)
2. **HIGH**: Fix generate_timeseries.R Rseq growth
3. **MEDIUM**: Fix stat_analysis.R list growth
4. **LOW**: Fix visualization list growth (only runs once)

---

## 📝 Additional Optimizations (Future)

### Nested Loop in generate_netdyn.R
Lines 62-68: Nested loops for computing Wcont could be vectorized using matrix operations.

### Nested Loop in estimate_MSAR.R
Lines 226-235: Triple nested loop for extracting/thresholding networks could be partially vectorized.

### Potential Parallelization
- Time series generation: Each of n_ts can be generated independently
- Model estimation: Each time series can be estimated independently
- Would require parallel package or future package

---

## ✅ Implementation Status

All critical performance optimizations have been **COMPLETED**:

### ✅ Fix #1: estimate_MSAR.R (COMPLETED)
**File**: `R/estimation/estimate_MSAR.R`
- Pre-allocated 10 vectors to maximum size (lines 153-166)
- Replaced all `c()` vector growth with direct indexing
- Multi-regime case optimized (lines 330-383)
- Single-regime case optimized (lines 405-454)
- Added vector trimming before summarize_cor() (lines 483-494)
- **Status**: ✅ Fully implemented

### ✅ Fix #2: generate_timeseries.R (COMPLETED)
**File**: `R/generation/generate_timeseries.R`
- Pre-allocated Rseq vector to known size (line 188)
- Replaced `Rseq <- c(Rseq, reg_index)` with `Rseq[m-1] <- reg_index` (line 202)
- **Status**: ✅ Fully implemented

### ✅ Fix #3: List Pre-allocation (COMPLETED)
**File**: `R/analysis/stat_analysis.R`
- Pre-allocated corr_list to estimated max size (lines 5-9)
- Replaced `append(corr_list, list(temp))` with direct indexing (line 31)
- Added list trimming after loop (line 40)
- **Status**: ✅ Fully implemented

**File**: `R/visualization/Create_3D_surfaceplots_for_means.R`
- Pre-allocated data_list to exact size (lines 5-12)
- Pre-allocated subplot_list in plotting loop (lines 210-212)
- Replaced all `append()` calls with direct indexing
- **Status**: ✅ Fully implemented

---

## 🔧 Next Steps

1. ✅ ~~Implement Fix #1 (estimate_MSAR.R)~~
2. ✅ ~~Implement Fix #2 (generate_timeseries.R)~~
3. ✅ ~~Implement Fix #3 (stat_analysis.R and visualization)~~
4. **TODO**: Run profiling script to benchmark improvements
5. **TODO**: Compare before/after timing results
6. **FUTURE**: Consider parallelization for production runs (using `parallel` or `future` package)
