# appa genotype effect on pE11 and 4G8 immunoreactivity

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23138739.svg)](https://doi.org/10.5281/zenodo.23138739)

Data and code for the Bayesian analysis of pE11 and 4G8 immunoreactivity in
6-month-old *appa*+/+ and *appa*−/− turquoise killifish brains (de Bakker et al., *Nature Aging*, 2026).

## Data

`data/<antibody>_intensity_<region>.csv`: mean immunofluorescence intensity in the
anterior rhombencephalon (AR) and hypothalamus superior lobe (HSL). Each row is one
slide carrying one *appa*+/+ section (column 1) and one *appa*−/− section (column 2),
stained and imaged together.

| Antibody | AR slides | HSL slides |
|---|---|---|
| pE11 | 6 | 5 |
| 4G8  | 6 | 5 |

## Model

`stan/appa.stan`, fitted separately for each antibody:

```
log(intensity) ~ Normal(mu[region] + alpha[slide] + delta[region] * WT, sigma)
mu ~ N(0, 2);  alpha ~ N(0, 1);  delta ~ N(0, 1);  sigma ~ Exponential(1)
```

`delta[region]` is the genotype effect as a log WT/KO ratio; `exp(delta)` is the fold
difference. For pE11, every slide has its own effect. For 4G8, slides are numbered
within each region, so slide *k* in AR and slide *k* in HSL share one slide effect,
as in the published fit.

## Run

From the repository root (R ≥ 4.3, rstan, tidyverse, ggdist):

```sh
Rscript R/fit_genotype_effect.R
```

This fits both models (4 chains × 2,000 iterations, 1,000 warm-up, `adapt_delta = 0.95`,
seed 1234) and writes to `results/`:

- `fit_pE11.rds`, `fit_4G8.rds`: stanfit objects
- `genotype_effect_summary.csv`: posterior median, 90% credible interval, P(WT > KO) and fold difference per region
- `genotype_effect_<antibody>.pdf`: posterior distributions of the genotype effect
- `sessionInfo.txt`: software versions

## Results

| Antibody | Region | log(WT/KO), median | 90% credible interval | P(WT > KO) |
|---|---|---|---|---|
| pE11 | AR  | 0.86 | [0.54, 1.18]  | >0.999 |
| pE11 | HSL | 0.76 | [0.41, 1.12]  | 0.998 |
| 4G8  | AR  | 0.34 | [−0.10, 0.78] | 0.90 |
| 4G8  | HSL | 0.47 | [−0.01, 0.94] | 0.95 |

All chains converged (R̂ ≤ 1.005, no divergent transitions).

## Citation

Archived on Zenodo: https://doi.org/10.5281/zenodo.23138739 (all versions; v1.0.0: https://doi.org/10.5281/zenodo.23138740).
