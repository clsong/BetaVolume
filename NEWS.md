# BetaVolume 0.2.0

This release aligns the package with the published method (Song et al. 2025,
*Ecological Monographs* 95:e70008, <https://doi.org/10.1002/ecm.70008>),
fixes several bugs, makes the package much faster and lighter, and adds full
documentation, a vignette and a test suite. A companion Julia package lives in
`julia/`.

## Changes that can alter results

* **New `weight_scheme` argument for the weighted embedding.** With
  `weights = TRUE`, a composition occurring `n_k` times among `N` communities
  with `m` unique compositions is weighted by either
  - `"legacy"` (default): `N * n_k / sum(n_l^2)`, the weights of version 0.1.0,
    used for all analyses of the paper (Figs. 9 and 10); or
  - `"paper"`: `m * n_k / N`, as printed in Eq. 5 and Fig. 4 of the paper; or
  - `"bounded"`: `n_k / max(n_l)`. No weight exceeds one, so weighting never
    increases beta diversity and `beta_vol <= d` is guaranteed for
    presence/absence data.

  The default is unchanged from version 0.1.0, so `betavolume(weights = TRUE)`
  returns the same values as before. The `"paper"` weights are larger by the
  factor `m * sum(n_l^2) / N^2 >= 1` and can push beta diversity above its
  nominal maximum `d` (about one in ten random presence/absence
  metacommunities; 2.12 for the two-species example of Fig. 4), which is why
  they are not the default. Contributions are identical under all schemes.
* **`disentangle_individual_contribution(weights = TRUE)` treats each unique
  composition as one unit**, as in Fig. 6 of the paper. Version 0.1.0 removed
  one row at a time, so every duplicated composition got a contribution of
  zero. The result is now named by the rows sharing each composition (for
  example `"2, 3"`).
* **`get_compressed_angle()` is normalised so that the identity matrix gives
  1** (the factor `2^d` of Eq. 8 was missing). Version 0.1.0 returned
  `Omega^(1/d) / 2`, which halved beta diversity when passed to `betavolume()`.

## Bug fixes

* The exact (`"deterministic"`) method no longer silently returns 0 on
  degenerate presence/absence data. Qhull's default `"Tv"` option aborts with
  an internal error on such inputs, and the error was caught and turned into
  0. Flat hulls are now detected from the rank of the point cloud and genuine
  Qhull errors are reported.
* `get_compressed_angle()` failed on R >= 4.2 with "the condition has length
  > 1".
* `get_compressed_angle()` used an absolute error tolerance (0.001) for an
  orthant probability of order `2^-d`, so it was inaccurate for ten or more
  species. It now uses closed forms for up to three species and a relative
  tolerance otherwise.
* `beta_temporal_overlap()` did not include the origin in the two hulls
  (Eq. 6), so it failed on the example of Fig. 5 of the paper.
* The hyper-normal approximation returned `NaN` when rounding produced tiny
  negative eigenvalues. Singular covariance matrices now give 0.
* No more deprecation warnings from `dplyr` on every weighted call.

## New features

* `beta_temporal_overlap(oriented = TRUE)` implements the oriented overlap of
  Fig. 5, which distinguishes synchronous from asynchronous changes.
* `disentangle_individual_contribution()` gains `hypervolume_method`,
  `dim_threshold` and `preprocess` arguments, and returns a named vector.
* `hypervolume_method = "exact"` is accepted as an alias of
  `"deterministic"`.
* Data frames, logical matrices and sparse matrices are accepted everywhere,
  and inputs are validated with informative error messages. Sparse matrices
  are detected automatically, so `sparsity` is no longer needed.

## Performance

* Contributions with the exact method recompute the hull only after removing
  hull vertices. With the hyper-normal method, all leave-one-out
  determinants come from a single Cholesky factorisation (matrix determinant
  lemma). Both are typically one to three orders of magnitude faster.
* Duplicated compositions are grouped with vectorised base R instead of
  `dplyr`, and duplicate points are removed before calling Qhull.
* The hyper-normal approximation uses a Cholesky factorisation instead of an
  eigendecomposition.

## Dependencies and documentation

* Dropped the dependencies on dplyr, tidyr, tibble, purrr and rlang. The
  package now imports only geometry, Matrix, mvtnorm, magrittr and base
  packages.
* All help pages were rewritten with the formulas of the paper, details,
  references and runnable examples. A new vignette,
  `vignette("BetaVolume")`, reproduces the worked examples of the paper.
* A testthat suite checks the package against the paper's worked examples,
  independent geometric computations, and results of version 0.1.0.
* Added `citation("BetaVolume")`.

# BetaVolume 0.1.0

* Initial version.
