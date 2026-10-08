# Bayesian Hierarchical Modelling in Stan: A Worked Progression

Fifteen Stan models, built up from a one-parameter beta-binomial to a negative binomial regression with correlated random intercepts and slopes, with an R workbook that fits each one through CmdStanR.

This is a self-study project I use to practise the full Bayesian workflow outside my research. Every experiment follows the same loop: simulate data from known parameters in a clinical-style scenario, fit the model, check the sampler, check the model against the data, and compare against the truth and a frequentist fit. All data are simulated; nothing here comes from a real study.

## What it covers

| # | Scenario | Models | Ideas practised |
|---|---|---|---|
| 1 | Response rate in a small trial | `beta_binom` | Conjugate priors, prior sensitivity, posterior probability of exceeding a threshold |
| 2 | Disease incidence with exposure | `poisson_gama`, `poisson_gamma_pooled` | Rates with person-time, pooling across counties of very different size |
| 3 | Blood pressure regression | `normal_regression`, `normal_regression_p2` | Comparison with OLS, centring predictors, posterior predictive bands |
| 4 | ICU mortality vs. severity score | `logistic_reg` | Near-separation, flat versus weakly informative priors |
| 5 | Over-dispersed counts | `poisson_reg`, `negbinom` | Posterior predictive checks on variance/mean, maximum, zeros; LOO comparison |
| 6 | Multi-site means | `no_pool`, `complete_pool`, `hierarchical` | Partial pooling, shrinkage factors, non-centred parameterisation |
| 7 | Hospital mortality profiling | `hier_logistic` | Random intercepts on the logit scale, intraclass correlation, funnel plots, raw versus adjusted rates |
| 10 | Readmission counts across hospitals | `m1_pooled_nb`, `m2_hier_int_nb`, `m3_hier_slopes_nb` | Prior predictive simulation, correlated random effects with an LKJ Cholesky prior, divergence and tree-depth checks, PPCs by group, LOO with Pareto-k diagnostics |

The numbering follows my study notes, so it is not contiguous.

The last model in the sequence, `m3_hier_slopes_nb.stan`, shows where the progression ends up:

```stan
parameters {
    real mu_alpha;
    real mu_beta;
    vector<lower=0>[2] tau;
    real<lower=0> phi;
    cholesky_factor_corr[2] L_Omega;
    matrix[2, J] z;
}
transformed parameters {
    matrix[2, J] u = diag_matrix(tau) * L_Omega * z;   // non-centred
    vector[J] alpha  = mu_alpha + u[1, ]';
    vector[J] beta_j = mu_beta  + u[2, ]';
}
model {
    L_Omega ~ lkj_corr_cholesky(2);
    to_vector(z) ~ normal(0, 1);
    y ~ neg_binomial_2_log(alpha[hospital] + beta_j[hospital] .* severity, phi);
    // priors on mu_alpha, mu_beta, tau, phi omitted here
}
```

Every model also writes `y_rep` and `log_lik` in `generated quantities`, so posterior predictive checks and PSIS-LOO work the same way throughout.

## Repository layout

```
bayesian_workflow.qmd   R/Quarto workbook: simulation, fitting, diagnostics and plots for every experiment
stan/                   the 15 Stan programs
```

## Running

Requires R, [CmdStan](https://mc-stan.org/users/interfaces/cmdstan) and the packages below.

```r
install.packages(c("tidyverse", "posterior", "bayesplot", "loo"))
install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
cmdstanr::install_cmdstan()
```

Open `bayesian_workflow.qmd` and run the chunks for an experiment in order; each one is self-contained and sets its own seed. Model paths are relative to the repository root.

## Skills demonstrated

Stan · CmdStanR · hierarchical and multilevel models · non-centred parameterisation · LKJ priors for correlated random effects · prior and posterior predictive checking (bayesplot) · model comparison with PSIS-LOO · MCMC diagnostics · simulation-based validation · R and ggplot2

## Author

Nishith Ranjon Roy (Pritom) · [LinkedIn](https://www.linkedin.com/in/pritom65/) · [Google Scholar](https://scholar.google.com.pk/citations?user=HjqFwTUAAAAJ&hl=en)
