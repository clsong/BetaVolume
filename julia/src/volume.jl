# Hypervolume estimators and betavolume().
#
# Internal functions take the *preprocessed* embedding (rows = points, columns = axes);
# the origin is always added here (it is part of the geometric object P, Eq. 3).

"""
Covariance matrix (denominator n - 1) of the embedded points plus the origin. For sparse
input it is formed as `(X'X - n μμ') / (n - 1)` so that the data are never densified.
"""
function _embedding_cov(X::Matrix{Float64})
    n0, d = size(X)
    n = n0 + 1
    μ = vec(sum(X; dims = 1)) ./ n
    Xc = X .- μ'
    C = Xc' * Xc
    C .+= μ * μ'           # the origin contributes (0 - μ)(0 - μ)'
    C ./= (n - 1)
    return (C .+ C') ./ 2
end

function _embedding_cov(X::SparseMatrixCSC{Float64,Int})
    n0, d = size(X)
    n = n0 + 1
    μ = vec(sum(X; dims = 1)) ./ n
    XtX = Matrix(X' * X)
    C = (XtX .- n .* (μ * μ')) ./ (n - 1)
    return (C .+ C') ./ 2
end

"""
Cholesky-based log-determinant of a symmetric PSD matrix; `-Inf` when the matrix is not
numerically positive definite (Cholesky fails, or a diagonal entry of the factor is
≤ sqrt(eps) × the largest one). Returns `(logdet, factorisation_or_nothing)`.
"""
function _logdet_psd(C::Matrix{Float64})
    isempty(C) && return (0.0, nothing)
    F = cholesky(Symmetric(C, :U); check = false)
    issuccess(F) || return (-Inf, nothing)
    dg = diag(F.U)
    mx = maximum(dg)
    any(<=(sqrt(eps()) * mx), dg) && return (-Inf, nothing)
    return (2 * sum(log, dg), F)
end

"Calibration constant of the hyper-normal estimator, 1 / Var(Bernoulli(1/2))."
const HYPER_NORMAL_CONSTANT = 4.0

function _hn_beta(X::DenseOrSparse)
    d = size(X, 2)
    ld, _ = _logdet_psd(_embedding_cov(X))
    isfinite(ld) || return 0.0
    return HYPER_NORMAL_CONSTANT * d * exp(ld / d)
end

function _exact_beta(X::DenseOrSparse)
    d = size(X, 2)
    return d * _hull_volume(Matrix(X))^(1 / d)
end

_beta_from_embedding(X, method::Symbol) = method === :exact ? _exact_beta(X) : _hn_beta(X)

function _check_angle(a)
    (a isa Real && isfinite(a) && a >= 0) ||
        throw(ArgumentError("`compressed_angle` must be a single finite non-negative number."))
    return Float64(a)
end

"""
    betavolume(Z; preprocess = true, weights = false, remove_unique = false,
               method = :auto, dim_threshold = 10, compressed_angle = 1.0,
               weight_scheme = :paper) -> Float64

Geometric beta diversity β_vol of Song et al. (2025). The metacommunity matrix is
embedded as a cloud of points (one per community or per species, whichever is more
numerous) in a `d`-dimensional space, `d = min(γ, N)`. Together with the origin the points
span a convex polytope `P` and

    β_vol = d ⋅ vol(P)^(1/d),

an effective number of communities between 0 (all communities identical) and `d` (for
presence/absence data: all 2^d − 1 non-empty compositions present). The origin is always
part of `P`, so empty communities (and species absent everywhere) do not change β.

# Methods (`method`)
- `:exact` (alias `:deterministic`): exact volume of the convex hull (Qhull, options
  `"Qt"`); exponential in `d`, practical for `d ≲ 15`. A hull that is not full dimensional
  has volume 0.
- `:hyper_normal`: ellipsoid approximation from the covariance `Σ` of the points plus the
  origin (denominator n − 1),

      β ≈ 4 d det(Σ)^(1/d),

  where 4 = 1/Var(Bernoulli(1/2)) calibrates the estimator so that the metacommunity with
  all 2^d compositions has β ≈ d. It scales like a variance (β(cZ) = c² β(Z)), whereas
  β_vol is linear (β(cZ) = c β(Z)). A singular `Σ` gives 0.
- `:auto` (default): `:exact` if `d ≤ dim_threshold`, `:hyper_normal` otherwise.

# Other arguments
- `preprocess`: pass `Z` through [`preprocess`](@ref) (orientation, `weights`,
  `remove_unique`, `weight_scheme`). With `false`, `Z` must already be an embedding
  (rows = points, columns = axes).
- `compressed_angle`: multiplicative factor Ω(S)^(1/d) accounting for species similarity
  or functional complementarity (Eq. 9), usually from [`compressed_angle`](@ref):
  `β_vol = d (Ω(S) vol(P))^(1/d)`.

Sparse matrices (`SparseMatrixCSC`) are supported; the hyper-normal covariance is then
computed without densifying the data.

# Examples
```julia
betavolume([1 0; 0 1])            # Fig. 1, metacommunity I:  2√(1/2) ≈ 1.414
betavolume([1 0; 0 1; 1 1])       # Fig. 1, metacommunity II: 2
fig4 = [1 0; 0 1; 0 1; 1 1; 1 1; 1 1]
betavolume(fig4; weights = true)  # 2√1.125
```
"""
function betavolume(Z; preprocess::Bool = true, weights::Bool = false,
                    remove_unique::Bool = false, method = :auto, dim_threshold::Real = 10,
                    compressed_angle::Real = 1.0, weight_scheme = :paper)
    angle = _check_angle(compressed_angle)
    scheme = _resolve_scheme(weight_scheme)
    X = preprocess ? _preprocess(Z, weights, remove_unique, scheme).points :
        _as_matrix(Z, "Z")
    n, d = size(X)
    m = _resolve_method(method, d, dim_threshold)
    n == 0 && return 0.0
    return angle * _beta_from_embedding(X, m)
end
