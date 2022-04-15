
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BetaVolume

<!-- badges: start -->
<!-- badges: end -->

The goal of BetaVolume is to …

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
#> ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.1 ──
#> ✓ ggplot2 3.3.5     ✓ purrr   0.3.4
#> ✓ tibble  3.1.6     ✓ dplyr   1.0.8
#> ✓ tidyr   1.2.0     ✓ stringr 1.4.0
#> ✓ readr   2.1.2     ✓ forcats 0.5.1
#> ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
#> x dplyr::filter() masks stats::filter()
#> x dplyr::lag()    masks stats::lag()
```

We analyze an empirical metacommunity collected collected by A. H.
Gentry and numerous additional collectors. This metacommunity has 81
species and 10 patches.

``` r
head(example_meta_composition)
#>                                        1 2 3 4 5 6 7 8 9 10
#> ANNONACEAE_ANOMIANTHUS_DULCIS          0 1 1 0 0 0 0 0 0  0
#> ANNONACEAE_DESMOS_M1                   0 1 0 0 0 0 0 0 0  0
#> ANNONACEAE_MITREPHORA_M1               0 0 0 0 0 0 0 0 0  1
#> APOCYNACEAE_AGANOSMA_MARGINATA         1 0 0 0 0 0 0 1 0  0
#> APOCYNACEAE_CHILOCARPUS_M1             1 0 0 0 0 0 0 0 0  0
#> APOCYNACEAE_HOLARRHENA_ANTIDYSENTERICA 0 0 0 0 0 0 0 0 1  0
```

As this data only contains information of presence/absence of species,
we adopt the duplication schemes to compute hypervolume beta diversity.
The hypervolume beta diversity of this metacommunity is

``` r
betavolume(example_meta_composition, weights = T)
#> [1] 3.099582
```

We can disentagle the individual contribution of a patch to the overall
beta diversity. We find that

``` r
disentangle_individual_contribution(example_meta_composition, weights = T) %>% 
  enframe(name = "species", value = 'contribution') %>% 
  arrange()
#> # A tibble: 81 × 2
#>    species contribution
#>      <int>        <dbl>
#>  1       1     5.26e-16
#>  2       2     1.18e-15
#>  3       3    -6.57e-16
#>  4       4     7.89e-16
#>  5       5     0       
#>  6       6     0       
#>  7       7     0       
#>  8       8     0       
#>  9       9     0       
#> 10      10     6.20e- 5
#> # … with 71 more rows
```
