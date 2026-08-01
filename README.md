# Bivariate censored extremes

R code accompanying *A functional limit theorem for heavy-tailed bivariate survival estimators* by Martin Bladt and Armelle Guillou.

The code implements tail inference for a heavy-tailed bivariate survival function under random right-censoring. It includes the Dabrowska product-limit estimator, pointwise Gaussian confidence intervals, adaptive threshold selection, simulation studies, and an insurance application.

## Files

- `functions.R`: estimation, inference, simulation, and plotting functions.
- `simulations.R`: Monte Carlo study and tail-set illustrations.
- `loss_alae.R`: application to the `lossalaefull` data from `CASdatasets`.

## Run

```sh
Rscript simulations.R
Rscript loss_alae.R
Rscript uk_data.R
```

The full simulation uses 1,000 replications in each of 24 scenarios and may take substantial time. Figures and LaTeX tables are written to `figures/` and `tables/`.
