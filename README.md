
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BetaVolume <img src="man/figures/logo_2.png" align="right" height="137" />

The goal of BetaVolume is to provides a new measure of beta diversity
based on the geometric embedding of metacommunities.

## Installation

You can install the development version of BetaVolume from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("clsong/betavolume")
```

## Example

Here we show a basic example of this package.

``` r
library(BetaVolume)
library(tidyverse)
```

We analyze an empirical metacommunity collected collected by A. H.
Gentry and numerous additional collectors. This metacommunity has 37
species and 10 patches.

``` r
head(example_meta_composition)
#>                                       1 2 3 4 5 6 7 8 9 10
#> ANACARDIACEAE_BUCHANANIA_RETICULATA   0 0 0 0 0 0 0 1 0  0
#> BURSERACEAE_CANARIUM_SUBULATUM        0 0 0 0 0 0 1 0 0  0
#> CLUSIACEAE_BAUHINIA_MALABARICA        0 0 0 0 1 0 0 0 0  0
#> CLUSIACEAE_CRATOXYLON_COCHINCHINENSIS 0 1 0 1 0 0 0 0 1  0
#> CLUSIACEAE_DALBERGIA_CULTRATA         0 0 0 1 0 0 0 0 0  0
#> CLUSIACEAE_DALBERGIA_NIGRESCENS       0 0 0 0 0 0 0 1 0  0
```

As this data only contains information of presence/absence of species,
we adopt the duplication schemes to compute hypervolume beta diversity.
The hypervolume beta diversity of this metacommunity is

``` r
betavolume(example_meta_composition, weights = T)
#> [1] 2.461931
```

We can disentangle the individual contribution of a patch to the overall
beta diversity. We find that

``` r
disentangle_individual_contribution(example_meta_composition, weights = T) %>% 
  enframe(name = "species", value = 'contribution') %>% 
  arrange(-contribution)
#> # A tibble: 37 × 2
#>    species contribution
#>      <int>        <dbl>
#>  1      10       0.211 
#>  2      21       0.140 
#>  3      29       0.122 
#>  4      22       0.111 
#>  5       9       0.0935
#>  6      12       0.0906
#>  7      30       0.0745
#>  8      17       0.0396
#>  9      27       0.0205
#> 10       4       0.0199
#> # … with 27 more rows
```
