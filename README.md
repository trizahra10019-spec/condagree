# condagree

**Ability-conditioned inter-rater agreement for constructed-response items**

Scoring reliability for constructed-response items is usually reported as one
number per item or test: percent exact agreement or Cohen's kappa between an
operational and a second rater. A single index cannot show *for whom* scoring
is unreliable. `condagree` estimates rater agreement as a function of
examinee ability, tests whether it is uniform, and separates two explanations
of an ability gradient:

* **score opportunity**: examinees at different ability levels produce
  different score categories, and partial-credit responses are intrinsically
  harder to score; and
* **residual ability effect**: scoring error changes with ability even for
  responses of the same correct category.

## Installation

```r
# install.packages("remotes")
remotes::install_github("trizahra10019-spec/condagree")
```

The package depends only on base R (`stats`, `graphics`, `grDevices`,
`splines`).

## Main functions

| Function | Purpose |
|---|---|
| `ca_data()` | Validate long-format double-scored data (one row per response) |
| `cond_agree()` | Kernel- or band-conditional exact agreement and item-stratified kappa, school-cluster bootstrap, Rubin pooling over plausible values |
| `agree_test()` | Test of uniform agreement: 1-df trend and spline omnibus test, item fixed effects, cluster-robust SEs, D1 pooling |
| `rater_model()` | Latent-category rater error model: score-opportunity intercepts (shared, by maximum score, or item-specific via `alpha_by`, compared by AIC), residual ability effect (gamma), directional severity, model-based decomposition |
| `severity_test()` | Direction of disagreement (second rater lower vs. higher) across ability |
| `simulate_agree()` | Data generator for simulation studies |

## Example

```r
library(condagree)
d <- simulate_agree(n_schools = 60, n_per_school = 20, gamma = 0.3,
                    alpha = c(zero = -4.5, partial = -2.2, full = -3.5), seed = 1)
x <- ca_data(d, first = "first", second = "second", item = "item",
             ability = paste0("PV", 1:5), cluster = "school", id = "student")

ca <- cond_agree(x, B = 200)   # conditional agreement and kappa
plot(ca)
agree_test(x)                  # is agreement uniform across ability?
rm <- rater_model(x, alpha_by = "maxscore")  # score opportunity vs. residual ability effect
rm
plot(rm)
severity_test(x)
```

## Replication

`inst/replication/` contains the scripts for the simulation study
(`sim_run.R`, `sim_rerun_rater_model.R`, `sim_rerun_tests.R`, `sim_summarize.R`), the PIRLS 2011 Indonesia illustration
(`pirls_application.R`, `pirls_sensitivity_alpha.R`, `supplement_figures.R`), a simulation of confounding by item composition (`sim_item_heterogeneity.R`, `sim_item_heterogeneity_summarize.R`), and an independent check of the
key PIRLS statistics (`verify_pirls_independent.R`). PIRLS data are not redistributed; they can be
downloaded from the IEA data repository.

## Citation

Fitrisia, A., Ningsih, T. Z., & Putri, S. K. (2026). condagree: An R package for
ability-conditioned inter-rater agreement in constructed-response scoring.
Manuscript submitted for publication.

## License

GPL (>= 3)
