# BetaVolume 0.2.0

This release aligns the package with the published method (Song et al. 2025,
*Ecological Monographs* 95:e70008, <https://doi.org/10.1002/ecm.70008>),
fixes several bugs, makes the package much faster and lighter, and adds full
documentation, a vignette and a test suite. A companion Julia package lives in
`julia/`.

## Changes that can alter results

* **Weighted embedding now follows Eq. 5 of the paper.** With
  `weights = TRUE`, a composition occurring `n_k` times among `N` communities
  with `m` unique compositions is weighted by `m * n_k / N` (Fig. 4 of the
  paper). Version 0.1.0 used `N * n_k / sum(n_l^2)`, which is proportional
  but smaller by the factor `m * sum(n_l^2) / N^2 >= 1`. Contributions are
  unaffected, but `betavolume(weights = TRUE)` values are larger than before.
  Use `weight_scheme = "legacy"` to reproduce results of version 0.1.0 (for
  example Fig. 9 of the paper).
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
