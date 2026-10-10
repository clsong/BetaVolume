# BetaVolume.jl

Julia port of the R package [BetaVolume](https://github.com/clsong/BetaVolume)
(version 0.2.0): geometric (hypervolume) measures of beta diversity from

> Song, C., Lu, M., Bennett, J. R., Gilbert, B., Fortin, M.-J., & Gonzalez, A. (2025).
> A geometric approach to beta diversity. *Ecological Monographs*, 95(1), e70008.
> [doi:10.1002/ecm.70008](https://doi.org/10.1002/ecm.70008)

A metacommunity matrix (communities × species) is embedded as a cloud of points in a
`d`-dimensional space, `d = min(#species, #communities)`. Together with the origin the
points span a convex polytope `P`, and beta diversity is

    β_vol = d ⋅ vol(P)^(1/d),

an effective number of communities ranging from 0 (all communities identical) to `d`
(all `2^d − 1` presence/absence compositions present). The R package is the reference
implementation; this port reproduces its results to ~1e-8 relative accuracy, which the
test suite checks against values computed by R (`test/reference.json`).

## Installation

```julia
using Pkg
Pkg.add(url = "https://github.com/clsong/BetaVolume", subdir = "julia")
```

Requires Julia ≥ 1.10. Dependencies: [DirectQhull](https://github.com/JuhaHeiskala/DirectQhull.jl)
(Qhull convex hulls) and [MvNormalCDF](https://github.com/PharmCat/MvNormalCDF.jl)
(orthant probabilities).

## Quick start: the worked examples of the paper

```julia
using BetaVolume

# Fig. 1: two species A and B
betavolume([1 0; 0 1])          # metacommunity I  {A}, {B}:        2√(1/2) ≈ 1.414
betavolume([1 0; 0 1; 1 1])     # metacommunity II {A}, {B}, {A,B}: 2 (the maximum)

# Fig. 4: duplicated presence/absence compositions, weighted embedding
fig4 = [1 0; 0 1; 0 1; 1 1; 1 1; 1 1]
betavolume(fig4; weights = true)                          # 1.818 (default :legacy weights)
betavolume(fig4; weights = true, weight_scheme = :paper)  # 2.121 (Eq. 5 as printed)
betavolume(fig4; weights = true, weight_scheme = :bounded)  # 1.414 (≤ d guaranteed)

# Fig. 6: contributions of the three unique compositions (Eq. 7)
c = contributions(fig4; weights = true)
c.values   # [0.1617, 0.3725, 0.4658]
c.units    # [[1], [2, 3], [4, 5, 6]]  (rows of fig4 forming each unit)
contributions([0.5 0; 0 1; 1.5 1.5]).values  # same, from the weighted matrix of Fig. 6A

# Fig. 5: temporal change as hypervolume overlap (Eq. 6)
t0 = [0 1; 1 1]; t1 = [0 1; 1 0]; t2 = [1 0; 0 1]
temporal_overlap(t0, t1)                    # √(1/2) ≈ 0.707
temporal_overlap(t1, t2)                    # 1
temporal_overlap(t1, t2; oriented = true)   # √2 (synchronous, opposite orientation)

# Fig. 7: species similarity shrinks, complementarity expands the hypervolume (Eqs. 8-9)
S = [1 0.5; 0.5 1]
compressed_angle(S)                                         # 0.640
betavolume([1 0; 0 1; 1 1]; compressed_angle = compressed_angle(S))
compressed_angle([1 -0.5; -0.5 1])                          # > 1

# Large metacommunities: fast hyper-normal approximation, sparse input supported
using SparseArrays
Z = sprand(5_000, 300, 0.05) .> 0
betavolume(Z)                       # d = 300 > dim_threshold → :hyper_normal
contributions(Z; method = :hyper_normal)
```

## API

| Julia | R | Description |
|---|---|---|
| `betavolume(Z; preprocess=true, weights=false, remove_unique=false, method=:auto, dim_threshold=10, compressed_angle=1.0, weight_scheme=:legacy)` | `betavolume()` | Geometric beta diversity. `method` ∈ `:auto`, `:exact` (alias `:deterministic`), `:hyper_normal`; `:auto` is exact for `d ≤ dim_threshold`. |
| `preprocess(Z; weights=false, remove_unique=false, weight_scheme=:legacy)` → `(points, ids, transposed)` | `preprocess_meta_composition()` | Orientation (rows = points), duplicate removal, weighted embedding. `ids[i]` is the row of `Z` (column if `transposed`) of embedded point `i`. |
| `contributions(Z; weights=false, remove_unique=false, method=:auto, dim_threshold=10, preprocess=true, weight_scheme=:legacy)` → `(values, units)` | `disentangle_individual_contribution()` | Normalised leave-one-out losses (Eq. 7). `units[k]` lists the `ids` forming unit `k` (R uses names such as `"2, 3"`). |
| `compressed_angle(S; m=25_000*size(S,1), rng=Random.default_rng())` | `get_compressed_angle()` | `Ω(S)^(1/d)`; exact for `d ≤ 3`, quasi-Monte Carlo (Genz) with `m` points otherwise. |
| `temporal_overlap(Z1, Z2; oriented=false, preprocess=true, weights=false, weight_scheme=:legacy)` | `beta_temporal_overlap()` | `(vol(P₁ ∩ P₂)/vol(P₁))^(1/d)`, optionally oriented. |

The R names are exported as aliases (`preprocess_meta_composition`,
`disentangle_individual_contribution`, `get_compressed_angle`, `beta_temporal_overlap`).
Inputs may be any `AbstractMatrix{<:Real}` (`Bool`, integer, `BitMatrix`,
`SparseMatrixCSC`, ...). R arguments map as follows: `hypervolume_method = "..."` →
`method = :...`, `weight_scheme = "..."` → `weight_scheme = :...` (strings are also
accepted); R's `sparsity` argument is unnecessary (sparse input is detected by type).

### Formulas

- **Exact**: `β = d vol(P)^(1/d)`, `P` = convex hull of the points and the origin (Qhull,
  options `"Qt"`, as in R); 0 when the points do not span `R^d`.
- **Hyper-normal**: `β ≈ 4 d det(Σ)^(1/d)`, `Σ` = covariance (n − 1 denominator) of the
  points plus the origin; 0 when `Σ` is singular. Scales quadratically, `β(cZ) = c²β(Z)`.
- **Weights** (Eq. 5): composition with multiplicity `n_k`, `N` non-empty rows, `m` unique
  compositions: `N n_k / Σ n_l²` (`:legacy`, the default, used for the paper's analyses)
  `m n_k / N` (`:paper`, Eq. 5 as printed) or `n_k / max n_l` (`:bounded`, guarantees
  `β_vol ≤ d` for presence/absence data).
- **Contributions** (Eq. 7):
  `C_i = (vol(P₀)^(1/d) − vol(P₋ᵢ)^(1/d)) / Σ_j (vol(P₀)^(1/d) − vol(P₋ⱼ)^(1/d))`.
- **Compressed angle** (Eq. 8): `Ω(S) = 2^d P(Z ≥ 0)`, `Z ~ N(0, (S'S)⁻¹)`.
- **Temporal overlap** (Eq. 6): `(vol(P_{t+1} ∩ P_t)/vol(P_t))^(1/d)`; oriented variant
  for simplices uses `|s₁| + |s₂|` when the signed volumes have opposite signs.

## Differences from the R package

- Results are returned as plain numbers / `NamedTuple`s instead of named vectors; `ids`
  and `units` play the role of R's row names. With `remove_unique = true`, `ids` refer
  to the original rows/columns (R numbers the deduplicated rows when the input has no
  names).
- `compressed_angle` for `d > 3` uses MvNormalCDF.jl with a fixed number of points `m`
  instead of mvtnorm's adaptive relative tolerance; pass a seeded `rng` for
  reproducibility. Values agree with R within the Monte Carlo error.
- Identical rows are detected by exact floating-point equality (R compares 15-significant
  digit strings).
- Leave-one-out arithmetic uses Float64 throughout (an integer overflow for very large
  groups of duplicated rows was found this way and is fixed in the R package too).
- The intersection of two hulls (temporal overlap) is computed with a self-contained
  simplex method (Chebyshev centre) and Qhull's halfspace-intersection duality, instead
  of `geometry::intersectn` + lpSolve. Results agree to ~1e-10.

## Performance notes

- The exact estimator is dominated by Qhull. Qhull's default merging (`Qx`, automatic in
  ≥ 5 dimensions) is expensive in high dimension: e.g. 255 points in 8-d produce ~2e5
  triangulated facets and take ~30 s in both R and Julia. The package calls Qhull directly
  (exit code checked, diagnostics silenced, only the volume read back), avoiding the
  per-facet conversion of `DirectQhull.ConvexHull`.
- Exact contributions recompute one hull per hull vertex that occurs once. This loop is
  multithreaded: start Julia with `JULIA_NUM_THREADS=k` (results are identical to the
  serial ones; Qhull is used through its reentrant library).
- Hyper-normal contributions need one Cholesky factorisation for all units (matrix
  determinant lemma); sparse input is never densified (`(X'X − nμμ')/(n − 1)`).
- Timings on the same matrices (best of 5; R 4.3.3, Julia 1.13):

  | task | Julia, 1 thread | Julia, 4 threads | R |
  |:--|--:|--:|--:|
  | `betavolume` exact, 200 × 8 binary, weights | 1.27 s | 1.07 s | 0.91 s |
  | `betavolume` hyper-normal, 2000 × 200 | 7.3 ms | 3.3 ms | 33 ms |
  | contributions exact, 200 × 8 binary, weights | 139 s | 44.8 s | 123 s |
  | contributions hyper-normal, 2000 × 100 | 4.6 ms | 6.7 ms | 20 ms |

  The hyper-normal path is 3–10× faster than R; single exact hulls are bound by Qhull
  and run at the same speed in both languages, but the multithreaded leave-one-out
  loop makes exact contributions about 3× faster.
- `benchmark/benchmarks.jl` times the main operations; `benchmark/benchmarks.R` runs the
  R package on the same matrices.

## Tests

```julia
using Pkg; Pkg.test("BetaVolume")
```

The suite covers the paper's worked examples (Figs. 1, 2, 4, 5, 6), known shapes,
invariances, brute-force checks of the fast leave-one-out algorithms, sparse input,
compressed angles against closed forms and Monte Carlo, the LP solver, and
cross-validation against the R implementation (`test/reference.json`, regenerated with
`Rscript julia/test/generate_reference.R`).
