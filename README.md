# fabPrediction

An **R** package for onstructing frequentist prediction regions using indirect information.

## References

Bersson and Hoff (2023). [Frequentist Prediction Sets for Species Abundance using Indirect Information](https://arxiv.org/pdf/2311.15860.pdf).

Bersson and Hoff (2022). [Optimal Conformal Prediction for Small Areas](https://arxiv.org/pdf/2204.08122.pdf).

## Installation

``` r
devtools::install_github("betsybersson/fabPrediction")  
```

## Usage

To load the package:

``` r
library(fabPrediction)  
```

The two main functions are 
- `predictionInterval`, which constructs prediction intervals for a continuous response. This function can be used to construct nonparametric FAB or distance-to-average conformal intervals, or parametric normal or Bayesian intervals.
- `predictionSet`, which constructs prediction sets for a categorical counts response. This function can be used to construct nonparametric FAB or direct sets, or a parametric Bayesian set.

Construction of basic FAB prediction regions are demonstrated below. Please see the vignette for full package capabilities, including empirical Bayes procedures to obtain estimates of prior hyperparameters based on auxiliary data.

### Continuous Response

We wlil demonstrate usage on a random normal sample of length 10.

``` r
y = rnorm(10)
y
```

A FAB prediction interval with 1-`alpha` coverage can be constructed for these data based on a prior parameters `mu` and `tau2` from a Normal-Normal working model:

``` r
y_PI = predictionInterval(y, method = "FAB",
      alpha = .05,
      mu = 0, tau2 = 1/2)
y_PI$bounds
```

and plotted:

``` r
plot(y_PI)
```

### Categorical Response

We wlil demonstrate usage on a random multinomial sample for 10 categories based on a heterogeneous prior concentration `gamma`.

``` r
gamma = c(10:1)
y = c(rmultinom(1,20,rdirichlet(gamma)))
y
```

A FAB prediction set with 1-`alpha` coverage can be constructed for these data based on an estimate of the prior parameter `gamma` from a Multinomial-Dirichlet working model:

``` r
y_PS = predictionSet(y, method = "FAB",
      alpha = .15,
      gamma = gamma)
y_PS$set
```

And this prediction set can be plotted:

``` r
plot(y_PS)
```

<> ### Examples

<> - vignette
<> - [Radon prediction in Bersson and Hoff (2022)](https://github.com/betsybersson/fab_sap/blob/main/radon_prediction.R) - need to update code to use the package...




