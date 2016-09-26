## mfa

`mfa` is an R package implementing Gibbs sampling for a Bayesian hierarchichal mixture of factor analysers for inference of bifurcations in single-cell data.

### Installation

```r
devtools::install_github("kieranrcampbell/mfa")
```

### Usage

For a cell-by-gene matrix of expression Y, MFA can be envoked via

```r
m <- mfa(Y, # gene expression
        iter = 2000, # MCMC iterations
        thin = 1, # MCMC samples to thin
        burn = iter / 2, # Number of MCMC samples to throw away
        b = 2, # Number of branches to model
        pc_initialise = 1, # Which principal component to initialise pseudotimes to
        collapse = FALSE, # Collapsed Gibbs sampling of branch assignments
        seed = 123L, # Random seed to set
        eta_tilde = mean(Y) # Hyperparameter eta tilde
        )
```

### Authors

Kieran Campbell & Christopher Yau

Wellcome Trust Centre for Human Genetics, University of Oxford

### Artwork

Upcoming...