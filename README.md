# transGFM: Transfer Learning for Generalized Factor Models to obtain 'blessing of dimension'.

[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![CRAN status](https://www.r-pkg.org/badges/version/transGFM)](https://CRAN.R-project.org/package=transGFM)
[![](https://cranlogs.r-pkg.org/badges/grand-total/transGFM?color=orange)](https://cran.r-project.org/package=transGFM)

## Overview

The `transGFM` package implements transfer learning methods for generalized factor models with support for continuous, count (Poisson), and binary data types. The package provides functions for single and multiple source transfer learning, source detection to identify positive and negative transfer sources, and information criteria for rank selection. The methods are particularly useful for high-dimensional data analysis where auxiliary information from related source datasets can improve estimation efficiency in the target domain.

The package is now available on [CRAN](https://cran.r-project.org/package=transGFM) under the name **transGFM**.

## Installation

### From CRAN (stable version)

```r
install.packages("transGFM")

# Load the package
library(transGFM)
```

### From GitHub (development version)

You can also install the development version of transGFM from GitHub:

```r
# Install devtools if you haven't already
install.packages("devtools")

# Install transGFM from GitHub
devtools::install_github("zjwangATsu/transGFM")

# Load the package
library(transGFM)
```

## Key Features

- **Single Source Transfer Learning**: Transfer knowledge from one source dataset to improve target estimation
- **Multiple Source Transfer Learning**: Leverage multiple related source datasets simultaneously
- **Source Detection**: Automatically identify positive and negative transfer sources
- **Flexible Data Types**: Support for continuous (Gaussian), count (Poisson), and binary data
- **Missing Data Handling**: Automatically handles missing values in source and target datasets
- **Information Criteria**: Built-in IC1 and IC2 criteria for selecting the number of factors
- **Cross-Validation**: Automatic parameter selection via cross-validation

## Main Functions

### Core Functions
- `transGFM()`: Single source transfer learning for generalized factor models
- `transGFM_multi()`: Multiple source transfer learning with AD (Aggregation then Debiasing) and DA (Debiasing then Aggregation) methods

### Source Detection Functions
- `source_potential()`: Identify potential sources based on rank comparison using IC criterion
- `source_detection()`: Detect positive/negative transfer sources using ratio-based method

### Utility Functions
- `ic_criterion()`: Information criterion (IC1/IC2) for selecting number of factors
- `relative_error()`: Calculate relative Frobenius norm error between matrices
- `identify()`: Factor decomposition using MLE (Maximum Likelihood Estimation)

## Supported Data Types

The package supports three types of data:

| Data Type | Distribution | Use Case |
|-----------|-------------|----------|
| `continuous` | Gaussian | Real-valued data |
| `count` | Poisson | Count data |
| `binary` | Binomial | Binary outcomes |

## Quick Start

### Load the Package

```r
library(transGFM)
set.seed(2025)
```

### Example 1: Single Source Transfer Learning (Count Data)

```r
# Generate Poisson data

# Source data (100 × 100 with 10% missing)
n1 <- 100; p1 <- 100; r <- 2
F_source <- matrix(runif(n1 * r, min = -2, max = 2), n1, r)
B_source <- matrix(runif(p1 * r, min = -2, max = 2), p1, r)
M_source <- F_source %*% t(B_source)
lambda_source <- exp(M_source)
X_source <- matrix(rpois(n1 * p1, as.vector(lambda_source)), n1, p1)

# Add 10% missing values to source
n_missing <- floor(n1 * p1 * 0.1)
missing_idx <- sample(n1 * p1, n_missing)
X_source[missing_idx] <- NA

# Target data (50 × 50, complete)
n0 <- 50; p0 <- 50
M_target_true <- M_source[1:n0, 1:p0]
lambda_target <- exp(M_target_true)
X_target <- matrix(rpois(n0 * p0, as.vector(lambda_target)), n0, p0)

# Run transGFM
result <- transGFM(
  source_data = X_source,
  target_data = X_target,
  r = 2,
  data_type = "count",
  lambda_seq = seq(0, 5, by = 1),
  K_cv = 3,
  verbose = FALSE
)

# Check results
print(paste("Optimal lambda:", result$optimal_lambda))
print(paste("Relative error:", relative_error(result$M_trans, M_target_true)))
```

### Example 2: Multiple Source Transfer Learning

```r
# Generate 3 source datasets
n1 <- 100; p1 <- 100; r <- 2
source_list <- list()

F_s <- matrix(runif(n1 * r, min = -2, max = 2), n1, r)
B_s <- matrix(runif(p1 * r, min = -2, max = 2), p1, r)
M_s <- F_s %*% t(B_s)

for (s in 1:3) {
  X_s <- matrix(rpois(n1 * p1, exp(M_s)), n1, p1)
  
  # Add missing values
  missing_rate <- 0.1 + (s - 1) * 0.02
  n_missing <- floor(n1 * p1 * missing_rate)
  missing_idx <- sample(n1 * p1, n_missing)
  X_s[missing_idx] <- NA
  
  source_list[[s]] <- X_s
}

# Target data
n0 <- 50; p0 <- 50
M_target_true <- M_s[1:n0, 1:p0]
X_target <- matrix(rpois(n0 * p0, exp(M_target_true)), n0, p0)

# Run transGFM_multi with AD method
result_AD <- transGFM_multi(
  source_data_list = source_list,
  target_data = X_target,
  r = 2,
  data_type = "count",
  method = "AD",
  lambda_seq = seq(0, 5, by = 1),
  K_cv = 3,
  verbose = FALSE
)

# Run transGFM_multi with DA method
result_DA <- transGFM_multi(
  source_data_list = source_list,
  target_data = X_target,
  r = 2,
  data_type = "count",
  method = "DA",
  verbose = FALSE
)

# Compare results
print(paste("AD method error:", relative_error(result_AD$M_trans, M_target_true)))
print(paste("DA method error:", relative_error(result_DA$M_trans, M_target_true)))
```

### Example 3: Source Detection (Rank-Based)

```r
# Generate 5 sources with different ranks
n1 <- 100; p1 <- 100
source_list <- list()

# Sources 1-2: rank 2 (same as target)
r_s <- 2
F_s <- matrix(runif(n1 * r_s, min = -2, max = 2), n1, r_s)
B_s <- matrix(runif(p1 * r_s, min = -2, max = 2), p1, r_s)
M_s <- F_s %*% t(B_s)

for (s in 1:2) {
  X_s <- matrix(rpois(n1 * p1, exp(M_s)), n1, p1)
  n_missing <- floor(n1 * p1 * 0.1)
  missing_idx <- sample(n1 * p1, n_missing)
  X_s[missing_idx] <- NA
  source_list[[s]] <- X_s
}

# Sources 3-5: rank 3 (different from target)
for (s in 3:5) {
  r_s_nega <- 3
  F_s_nega <- matrix(runif(n1 * r_s_nega, min = -2, max = 2), n1, r_s_nega)
  B_s_nega <- matrix(runif(p1 * r_s_nega, min = -2, max = 2), p1, r_s_nega)
  M_s_nega <- F_s_nega %*% t(B_s_nega)
  X_s_nega <- matrix(rpois(n1 * p1, exp(M_s_nega)), n1, p1)
  
  n_missing <- floor(n1 * p1 * 0.1)
  missing_idx <- sample(n1 * p1, n_missing)
  X_s_nega[missing_idx] <- NA
  
  source_list[[s]] <- X_s_nega
}

# Target data: rank 2
n0 <- 50; p0 <- 50
M_target <- M_s[1:n0, 1:p0]
X_target <- matrix(rpois(n0 * p0, exp(M_target)), n0, p0)

# Identify potential sources
result <- source_potential(
  X_sources = source_list,
  X0 = X_target,
  r_max = 5,
  ic_type = "IC1",
  data_type = "count",
  verbose = TRUE
)

print(result$positive_potential_sources)  # Should be c(1, 2)
print(result$negative_sources)            # Should be c(3, 4, 5)
print(result$r_target)                    # Should be 2
print(result$r_sources)                   # Should be c(2, 2, 3, 3, 3)
```

### Example 4: Selecting Number of Factors

```r
# Generate Poisson data with known rank
n <- 50; p <- 50; r_true <- 2

# Generate true factors
F_true <- matrix(runif(n * r_true, min = -1, max = 1), n, r_true)
B_true <- matrix(runif(p * r_true, min = -1, max = 1), p, r_true)
M_true <- F_true %*% t(B_true)

# Generate Poisson observations
lambda <- exp(M_true)
X <- matrix(rpois(n * p, as.vector(lambda)), n, p)

# Add 10% missing values
n_missing <- floor(n * p * 0.1)
missing_idx <- sample(n * p, n_missing)
X[missing_idx] <- NA

# Use IC1 to select rank
result_IC1 <- ic_criterion(
  X = X,
  r_max = 5,
  ic_type = "IC1",
  data_type = "count",
  verbose = FALSE
)

print(paste("True rank:", r_true))
print(paste("Estimated rank (IC1):", result_IC1$r_hat))

# Use IC2 to select rank
result_IC2 <- ic_criterion(
  X = X,
  r_max = 5,
  ic_type = "IC2",
  data_type = "count",
  verbose = FALSE
)

print(paste("Estimated rank (IC2):", result_IC2$r_hat))
```

## Main Parameters

### Data and Model Specification

- **`source_data`** / **`source_data_list`**: Source dataset(s) (may contain NA for missing values)
- **`target_data`**: Target dataset (n0 × p0 matrix)
- **`r`**: Number of factors
- **`data_type`**: Type of data - `"continuous"`, `"count"`, or `"binary"`

### Transfer Learning Parameters

- **`lambda_seq`**: Sequence of lambda values for cross-validation
- **`K_cv`**: Number of folds for cross-validation (default: 5)
- **`method`**: Fusion method for multiple sources - `"AD"` (Aggregation then Debiasing) or `"DA"` (Debiasing then Aggregation)

### Algorithm Parameters

- **`C`**: CJMLE projection constant (if NULL, auto-calculated)
- **`max_iter_cjmle`**: Maximum iterations for CJMLE (default: 30)
- **`max_iter_refine`**: Maximum iterations for refinement (default: 30)
- **`max_iter_nuclear`**: Maximum iterations for nuclear norm optimization (default: 100)
- **`verbose`**: Print progress information (default: FALSE)

### Information Criterion Parameters

- **`r_max`**: Maximum number of factors to consider
- **`ic_type`**: IC criterion type - `"IC1"` or `"IC2"`

## Methodological Details

The transGFM package employs:

1. **CJMLE Algorithm**: Maximum Likelihood Estimation for factor models with missing data
2. **Transfer Learning Framework**: 
   - Single source: Direct transfer with cross-validated regularization
   - Multiple sources: AD (aggregate then debias) or DA (debias then aggregate) methods
3. **Source Detection**:
   - Rank-based: Compare estimated ranks using IC criterion
   - Ratio-based: Use likelihood ratio to identify beneficial sources
4. **Missing Data**: Automatically handles missing values through indicator matrix approach
5. **Convergence**: Monitors likelihood changes until convergence

## Model Characteristics

The transGFM package is specifically designed to handle:

1. **High-dimensional data** where the number of variables may exceed observations
2. **Missing data** in both source and target datasets
3. **Multiple data types** including continuous, count, and binary
4. **Heterogeneous sources** with different sample sizes and missing rates
5. **Large Feature** through bless of dimension

The transfer learning framework can significantly improve estimation accuracy when:
- Source and target domains share similar latent factor structures
- Source datasets have more observations or lower missing rates
- Multiple complementary sources are available

## Dependencies

- R (≥ 3.5.0)
- stats

## Citation

The related papers are currently under review. We will update the citation information once they are published.

## Author

**Zhijing Wang**  
Shanghai Jiao Tong University (SJTU)  
Email: wangzhijing@sjtu.edu.cn

## Bug Reports and Issues

Please report any bugs or issues on the [GitHub Issues page](https://github.com/zjwangATsu/transGFM/issues). When reporting, include:

1. A clear description of the issue
2. Reproducible code example
3. Your session information (`sessionInfo()`)
4. Any relevant error messages

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request. For major changes, please open an issue first to discuss proposed modifications.

## License

This package is licensed under GPL-3. See the [LICENSE](LICENSE) file for details.
