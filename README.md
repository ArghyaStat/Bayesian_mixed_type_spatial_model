# Bayesian_mixed_type_spatial_model

This repository accompanies the article:

**“Scalable Bayesian inference for high-dimensional mixed-type multivariate spatial data”**

It contains all code and outputs required to reproduce the **simulation studies**, **real data analysis**, and **computational quality assessments** presented in the paper.

---

## Repository Structure

The repository is organized into three main components:

- **Simulated_data_analysis**
- **Real_data_analysis**
- **Quality_assesment**

Each folder directly corresponds to specific sections of the manuscript and is structured to facilitate reproducibility and clarity.

---

## 1. Simulated Data Analysis

📁 `Simulated_data_analysis/`

This folder contains the complete implementation of the simulation studies presented in the paper.

### Models Considered

The following mixed-type response settings are analyzed:

- **Binomial–Gaussian (BG)**
- **Binomial–Poisson (BP)**
- **Gaussian–Poisson (GP)**
- **Binomial–Gaussian–Poisson (BGP)**

### Experimental Design

For each response type, the simulations include:

- Sample sizes:  
  - $n = 100$ (small-scale)  
  - $n = 2500$ (large-scale)

- Spatial dependence parameter:  
  - $\phi = 0.1$ (weak dependence)  
  - $\phi = 0.3$ (strong dependence)

- Cross-response covariance structures:  
  - **Dependent $\bm{\Sigma}$**  
  - **Independent $\bm{\Sigma}$**

- Number of replications:  
  - **50 independent datasets**

### Model Comparisons

Each setting includes:

- **Joint model (proposed method)**
- **Separate model (baseline approach)**

### Contents

Each subfolder contains:

- Parallelized R scripts for simulation and inference
- MCMC implementation (warm-up and main sampling)
- Estimation summaries
- Predictive evaluation outputs

### Correspondence to Paper

- Estimation results → Tables (Simulation section)
- Predictive performance → ELPD / coverage tables
- Additional summaries → Appendix materials

---

## 2. Real Data Analysis

📁 `Real_data_analysis/`

This folder contains the real data application described in the paper.

### Subfolders

- `M1/`  
  Implementation of **Model 1** (baseline / structured specification)

- `M2/`  
  Implementation of **Model 2** (proposed joint framework)

- `Output_analysis/`  
  Post-processing of MCMC output, including:
  - Parameter summaries
  - Predictive performance
  - Tables and figures reported in the manuscript

### Correspondence to Paper

- Main results → Table 6 (and related figures)
- Supplementary outputs → Appendix

---

## 3. Quality Assessment

📁 `Quality_assesment/`

This folder contains additional experiments evaluating **computational efficiency, scalability, and sampler performance**, as discussed in the paper.

### Subfolders

#### (a) Vecchia Approximation Analysis

- `Vecchia_comparison/`  

Study of different conditioning sizes ($m$), including:

- Accuracy vs. computational cost trade-offs
- Impact on posterior inference

---

#### (b) Sampler Comparisons

Evaluation of different MCMC strategies:

- `joint_elliptical/` → Proposed sampler  
- `comp_wise_elliptical/` → Component-wise elliptical sampler  
- `compwise_pcn/` → Preconditioned Crank–Nicolson (pCN)  
- `joint_rwmh/` → Random Walk Metropolis–Hastings  

Each folder contains:

- Implementation of the sampler
- Performance metrics:
  - Effective sample size (ESS)
  - ESS per unit time
  - Mixing diagnostics

---

#### (c) Output Analysis

- Aggregated comparison results across samplers
- Figures and tables used in:
  - Sampler comparison sections
  - Appendix diagnostics

---

## Reproducibility Notes

- All simulations are **fully reproducible** using the provided scripts.
- Parallel computation is implemented using base R parallelization.
- No external libraries beyond base R (and `stats`) are required.

---

## Guidance for Users

To reproduce specific components:

- **Simulation results:**

- **Real data analysis:**

- **Methodological comparisons:**


---

## Citation

If you use this code, please cite:

*Scalable Bayesian inference for high-dimensional mixed-type multivariate spatial data*
