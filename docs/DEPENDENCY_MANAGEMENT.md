# Dependency Management Guide

## Overview

SimMSAR uses a centralized dependency management system to ensure all required R packages are available and loaded consistently across the project.

## Files

### R/dependencies.R
The main dependency management file that:
- Lists all required packages organized by category
- Provides functions to check, install, and load packages
- Automatically loads all packages when sourced

### DESCRIPTION
Standard R package metadata file that documents all dependencies for tools like `renv` or `packrat`.

## Usage

### In Your Scripts

Simply source the dependencies file at the beginning:

```r
source("R/dependencies.R")
```

This will automatically load all required packages.

### Check Package Status

```r
source("R/dependencies.R")
print_package_summary()
```

Output:
```
=== SimMSAR Package Dependencies ===

Core:                (3 packages)
  ✓ NHMSAR            (v1.0.0)
  ✓ netcontrol        (v0.2.0)
  ✓ huge              (v1.3.5)

Data:                (2 packages)
  ✓ dplyr             (v1.1.0)
  ✓ Matrix            (v1.5.1)
...
```

### Install Missing Packages

```r
source("R/dependencies.R")
install_missing_packages()
```

This will check for and install any missing packages.

### Install Specific Packages

```r
source("R/dependencies.R")
install_missing_packages(c("ggplot2", "dplyr"))
```

## Package Categories

### Core
Essential packages for MSAR functionality:
- **NHMSAR**: Non-homogeneous Markov-switching autoregressive models
- **netcontrol**: Network control analysis
- **huge**: High-dimensional graph estimation (nonparanormal transformation)

### Data Manipulation
- **dplyr**: Data manipulation and transformation
- **Matrix**: Sparse and dense matrix classes

### Statistical Modeling
- **lme4**: Linear mixed-effects models
- **lmerTest**: Tests for mixed models
- **lmPerm**: Permutation tests
- **dunn.test**: Dunn's test for multiple comparisons
- **effects**: Effect displays
- **sjPlot**: Statistical plots and tables

### Simulation
- **graphicalVAR**: Graphical vector autoregression
- **mvtnorm**: Multivariate normal distributions
- **abind**: Combine multi-dimensional arrays

### Regularization
- **lars**: LASSO and least angle regression
- **prettyGraphs**: Graph visualization

### Visualization
- **ggplot2**: Grammar of graphics plotting
- **plotly**: Interactive web-based graphs
- **cowplot**: Streamlined plot themes
- **RColorBrewer**: Color palettes

### Output/Reporting
- **knitr**: Dynamic report generation
- **kableExtra**: Complex table construction
- **xtable**: LaTeX/HTML table export

### Utilities
- **progress**: Progress bars for loops

## Best Practices

### DO:
✅ Source `R/dependencies.R` at the start of main scripts
✅ Document required packages in individual function files as comments
✅ Keep the dependency list in `R/dependencies.R` up to date
✅ Use `print_package_summary()` to verify all packages are installed

### DON'T:
❌ Add `library()` calls in individual function files
❌ Load packages multiple times
❌ Use `require()` instead of the centralized system
❌ Install packages in function files

## Migrating Existing Code

When adding new functions to the project:

1. **Remove library() calls** from individual files:
   ```r
   # OLD - Don't do this
   library(dplyr)
   library(ggplot2)

   # NEW - Add a comment instead
   # Dependencies are loaded centrally via R/dependencies.R
   # Required packages: dplyr, ggplot2
   ```

2. **Add new dependencies** to `R/dependencies.R`:
   ```r
   PACKAGES <- list(
     # ... existing packages ...

     new_category = c(
       "newpackage1",
       "newpackage2"
     )
   )
   ```

3. **Update DESCRIPTION** file:
   ```
   Imports:
       ...existing packages...,
       newpackage1,
       newpackage2
   ```

## Using with Package Managers

### renv (Recommended)

```r
# Initialize renv for the project
renv::init()

# Install all dependencies
source("R/dependencies.R")
install_missing_packages()

# Take a snapshot
renv::snapshot()

# Restore packages (e.g., on a new machine)
renv::restore()
```

### packrat

```r
# Initialize packrat
packrat::init()

# Install dependencies
source("R/dependencies.R")
install_missing_packages()

# Take snapshot
packrat::snapshot()
```

## Troubleshooting

### Package Won't Load

```r
# Check if package is installed
"packagename" %in% installed.packages()[, "Package"]

# Try manual installation
install.packages("packagename")

# Check for conflicts
conflicts(detail = TRUE)
```

### Version Conflicts

Some packages may require specific versions. Document these in DESCRIPTION:

```
Imports:
    packagename (>= 1.0.0)
```

### NHMSAR or netcontrol Not Found

These packages may not be on CRAN. Install from GitHub:

```r
# Install devtools if needed
install.packages("devtools")

# Install from GitHub (update with correct repository)
devtools::install_github("username/NHMSAR")
devtools::install_github("username/netcontrol")
```

## Future Enhancements

- [ ] Add version pinning for reproducibility
- [ ] Integrate with Docker for containerized environments
- [ ] Create pre-commit hooks to check dependencies
- [ ] Add CI/CD integration for automated dependency checks
