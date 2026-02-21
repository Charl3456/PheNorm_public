# EHR Phenotyping Simulation Study

Code repository for the simulation study described in:

> **Hong Y, Nelson JC, Williamson BD.** *Performance of weakly-supervised electronic health record-based phenotyping methods in rare-outcome settings.*

This repository reproduces all simulation results comparing weakly-supervised EHR phenotyping algorithms — **PheNorm**, **MAP**, and **sureLDA** — across data-generating mechanisms, outcome prevalence levels, and silver label informativeness, with a focus on rare-event settings relevant to vaccine safety surveillance.

---

## Background

PheNorm, MAP (Multimodal Automated Phenotyping), and sureLDA (Surrogate-guided Latent Dirichlet Allocation) are computable phenotyping methods that combine structured EHR features (ICD codes) with features derived from clinical notes via natural language processing (NLP). They are *weakly-supervised*, using silver-standard proxy labels rather than gold-standard chart-reviewed outcomes.

These methods were originally developed and evaluated for chronic, common conditions. This study assesses how they perform in rare-outcome settings such as vaccine adverse events (e.g., anaphylaxis), where outcome rates may be ~5%, documentation is inconsistent, and data sparsity is a key challenge.

---

## Repository Structure

| File | Description |
|---|---|
| `1_Helper_data_gen.R` | Data generation functions for all three data-generating mechanisms (simplified, sureLDA/LDA, complex) and preprocessing for MAP input |
| `2_Helper_data_process.R` | Preprocessing functions that format datasets for PheNorm, sureLDA, and sureLDA-MAP inputs |
| `3_Helper_run.R` | Model-fitting, prediction, and data-splitting functions for all algorithm variants |
| `4_run_simulation.R` | Main simulation runner — loops over replications, runs all algorithms, computes performance metrics, and saves results |
| `Combine_simulation_results.Rmd` | Aggregates per-run `.RData` files across parallel runs into a single combined dataset |
| `Table_Figure_generation.Rmd` | Produces all manuscript tables and figures |

---

## Algorithms

Eight algorithm versions are compared, varying initialization and tuning parameters:

| Algorithm | Description |
|---|---|
| **PheNorm-v1** | Standard PheNorm with fixed variance initialization (= 1) |
| **PheNorm-v2** | PheNorm with data-adaptive variance initialization (= half SD of normalized scores), as used internally by sureLDA |
| **MAP-v1** | Standard MAP fit on all patients with non-zero ICD codes |
| **MAP-v2** | MAP fit only on filter-positive patients, as implemented within sureLDA |
| **sureLDA-v1** | sureLDA with PheNorm-v1 priors |
| **sureLDA-v2** | sureLDA with MAP-v1 priors |
| **sureLDA-v3** | sureLDA with PheNorm-v2 priors |
| **sureLDA-v4** | sureLDA with MAP-v2 priors |
| **ICD Logit** | Logistic regression on ICD codes only (baseline) |

---

## Simulation Design

### Scenarios

Performance is evaluated under a **2×2 factorial design**:

| Dimension | Levels |
|---|---|
| **Outcome prevalence** | Rare (5%) vs. Common (40%) |
| **Silver label informativeness** | Highly predictive vs. Weakly predictive |

The 5% rare setting reflects vaccine adverse event rates; the 40% common setting reflects the chronic disease conditions for which these methods were originally designed.

### Data-Generating Mechanisms

Three mechanisms of increasing complexity are implemented:

**Simplified** — an idealized baseline where disease probabilities follow normal distributions and silver labels (ICD and NLP counts) follow Poisson mixtures with clear case-control separation (informative) or identical distributions (non-informative). Each patient contributes exactly one clinical note.

**sureLDA (LDA)** — implements the same generative model used in the original sureLDA paper. Latent continuous outcomes are thresholded to binary disease status; silver labels follow Gamma distributions scaled by healthcare utilization. An additional 150 auxiliary NLP variables (noise features) are included to simulate high-dimensional EHR data.

**Complex** — the most realistic mechanism. Patient demographics (age, sex, race, exposure history) influence disease probability through logistic regression. Silver label generation follows a realistic temporal sequence: text mentions → NLP features → ICD codes, creating feature dependencies characteristic of real-world EHR phenotyping.

### Evaluation

For each of the 3 × 4 = 12 data-generating mechanism × scenario combinations, **2,500 independent datasets** (n = 10,000 each) are generated. Each dataset is split into:
- A large **training set** (n = 9,800) with gold-standard labels hidden from algorithms
- A small labeled **test set** (n = 200)

Performance metrics computed on both sets:
- AUC (primary discrimination metric)
- F1 score, precision, recall, accuracy
- Probability MSE, MAE, and correlation vs. true latent probabilities (calibration)
- Optimal classification threshold (Youden index)

Results are reported as means and standard deviations over the 2,500 replications.

---

## How to Run

### 1. Install Required R Packages

```r
install.packages(c("dplyr", "Matrix", "pROC", "logger", "ggplot2",
                   "tidyr", "knitr", "kableExtra", "viridis", "purrr", "here"))
```

Install the phenotyping method packages (see each package's documentation for installation):

```r
# PheNorm, MAP, sureLDA
# (add GitHub/CRAN installation instructions here)
```

### 2. Run a Single Simulation Batch

Each call to `4_run_simulation.R` takes a run number as a command-line argument and executes 32 replications, saving results to `simulation_results7/`.

```bash
Rscript 4_run_simulation.R 1   # Run 1 — seed: 20250607
Rscript 4_run_simulation.R 2   # Run 2 — seed: 20251607
# ... repeat for as many runs as desired
```

Each run produces:
- `simulation_results7/simulation_results_narm_run<N>.RData`
- `simulation_results7/simulation_results_epsilon_run<N>.RData`
- `simulation_results7/diagnostics_run<N>.txt`

### 3. Run in Parallel

Runs are embarrassingly parallel with unique seeds (`20250606 + run_number × 1000`). On a cluster or multi-core machine:

```bash
for i in $(seq 1 50); do
  Rscript 4_run_simulation.R $i &
done
```

The study used ~50 runs × 32 replications = **~2,500 total replications** per scenario.

### 4. Handling Missing Predictions

Rare numerical instabilities in Bayesian updates (typically <5% of individuals, occurring when true probabilities approach zero) can produce undefined predictions. Two strategies are implemented in parallel:

- **`na.rm`** — imputes missing predictions with the column mean
- **`epsilon`** — replaces missing values with ε = 1×10⁻¹⁰ (primary strategy reported in the paper)

Both are saved separately; results were robust to the choice of strategy.

### 5. Combine Results

Update the directory paths at the top of `Combine_simulation_results.Rmd` to point to your local results folders, then knit:

```r
rmarkdown::render("Combine_simulation_results.Rmd")
```

### 6. Generate Tables and Figures

```r
rmarkdown::render("Table_Figure_generation.Rmd")
```

---

## Output Structure

Each `.RData` file contains:

- **`results_narm`** / **`results_epsilon`** — long-format data frame with one row per (replication × data type × scenario × algorithm × train/test set)
- **`summary_stats_narm`** / **`summary_stats_epsilon`** — grouped summary with mean and SD of all metrics across replications

---

## Key Findings (Summary)

- No single method dominated across all scenarios.
- Under **simplified data**, PheNorm and MAP achieved near-perfect AUC; sureLDA was unstable and prone to overfitting.
- Under **complex data with rare outcomes and informative labels** — the most practically relevant setting — sureLDA-v3 (sureLDA with PheNorm-v2 priors) achieved the highest AUC (0.95) and F1 (0.64); PheNorm-v2 had the best probability calibration (MSE = 0.14).
- Performance varied by up to 0.28 AUC points depending on tuning parameter choices, underscoring the importance of evaluating multiple configurations rather than relying on defaults.
- Algorithm-guided probability sampling for chart review selected encounters with higher rates of anaphylaxis-related NLP mentions and emergency treatments compared to simple random sampling.

---

## Funding

Supported by Task Order 75D30122F00001 from the U.S. Centers for Disease Control and Prevention (CDC). The content is solely the responsibility of the authors and does not necessarily represent the official views of the CDC or the U.S. Government.
