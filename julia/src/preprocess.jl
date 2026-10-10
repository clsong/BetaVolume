"""
    preprocess(Z; weights = false, remove_unique = false, weight_scheme = :legacy)
        -> (points, ids, transposed)

Prepare a metacommunity matrix for the geometric embedding of Song et al. (2025)
(R: `preprocess_meta_composition`):

1. optionally remove duplicated rows and then duplicated columns (`remove_unique = true`;
   despite the historical name, one copy of each unique row/column is kept);
2. orient the matrix so that **rows are the embedded points and columns the axes**: the
   matrix is transposed when it has fewer rows than columns, so that the embedding
   dimension is `d = min(#species, #communities)`;
3. optionally apply the weighted embedding of duplicated compositions (`weights = true`,
   Eq. 5 and Fig. 4 of the paper).

# Weighted embedding
All-zero rows are dropped and each row is multiplied by a weight proportional to the
multiplicity `n_k` of its composition. With `N` non-empty rows and `m` unique
compositions:

- `weight_scheme = :legacy` (default): `w_k = N n_k / Σ_l n_l²`. The weights average to
  one over communities. This is the normalisation used for all analyses of Song et al.
  (2025), including Figs. 9 and 10, and in practice it keeps `β_vol` within `[0, d]`.
- `weight_scheme = :paper`: `w_k = m n_k / N`, as printed in Eq. 5 and Fig. 4 of the
  paper (weights average to one over unique compositions). These weights can exceed
  one substantially, so `β_vol` can exceed its nominal maximum `d` (2.12 for the
  two-species example of Fig. 4).
- `weight_scheme = :bounded`: `w_k = n_k / max_l n_l`. No weight exceeds one, so the
  weighted hull lies inside the unweighted one: weighting never increases beta
  diversity and `β_vol ≤ d` is guaranteed for presence/absence data.

All three leave a metacommunity without duplicates unchanged and give identical relative
contributions ([`contributions`](@ref)), but different `betavolume(...; weights = true)`
values.

# Arguments
- `Z`: an `AbstractMatrix{<:Real}` (dense, `Bool`, integer or `SparseMatrixCSC`) with
  communities in rows and species in columns, or the reverse. Entries are presence/absence
  or any non-negative measure of importance (abundance, biomass, ...).

# Returns
A `NamedTuple`:
- `points`: `Matrix{Float64}` (or `SparseMatrixCSC{Float64,Int}` for sparse input) whose
  rows are the embedded points and whose columns are the axes (`nrow ≥ ncol`, except when
  `weights = true` dropped empty rows);
- `ids`: for each row of `points`, its index among the rows of the (possibly transposed)
  input `Z`, i.e. the row of `Z`, or the column of `Z` when `transposed`. This plays the
  role of the R row names (with `remove_unique = true` it refers to the *original*
  rows/columns);
- `transposed`: whether the input was transposed.

# Example
```julia
fig4 = [1 0; 0 1; 0 1; 1 1; 1 1; 1 1]
preprocess(fig4; weights = true).points   # rows scaled by 0.5, 1, 1, 1.5, 1.5, 1.5
```
"""
function preprocess(Z; weights::Bool = false, remove_unique::Bool = false,
                    weight_scheme = :legacy)
    return _preprocess(Z, weights, remove_unique, _resolve_scheme(weight_scheme))
end

"""
    preprocess_meta_composition(Z; weights = false, remove_unique = false, weight_scheme = :legacy)

Alias of [`preprocess`](@ref), with the name of the R function.
"""
const preprocess_meta_composition = preprocess

function _preprocess(Z, weights::Bool, remove_unique::Bool, scheme::Symbol,
                     arg::AbstractString = "Z")
    X = _as_matrix(Z, arg)
    return _preprocess_matrix(X, weights, remove_unique, scheme)
end

function _preprocess_matrix(X::M, weights::Bool, remove_unique::Bool,
                            scheme::Symbol) where {M<:DenseOrSparse}
    rows = collect(1:size(X, 1))
    cols = collect(1:size(X, 2))
    if remove_unique
        rows = _row_groups(X).first
        X = X[rows, :]
        keepc = _col_groups(X).first
        cols = cols[keepc]
        X = X[:, keepc]
    end
    transposed = size(X, 1) < size(X, 2)
    if transposed
        X = permutedims(X)
        ids = cols
    else
        ids = rows
    end
    if weights
        X, keep = _weight_composition(X, scheme)
        ids = ids[keep]
    end
    return (points = X::M, ids = ids, transposed = transposed)
end

function _nonzero_rows(X::Matrix{Float64})
    n, d = size(X)
    nz = falses(n)
    @inbounds for j in 1:d, i in 1:n
        X[i, j] != 0 && (nz[i] = true)
    end
    return nz
end

function _nonzero_rows(X::SparseMatrixCSC{Float64,Int})
    nz = falses(size(X, 1))
    rv = rowvals(X)
    vals = nonzeros(X)
    @inbounds for p in eachindex(vals)
        vals[p] != 0 && (nz[rv[p]] = true)
    end
    return nz
end

_scale_rows(X::Matrix{Float64}, w::Vector{Float64}) = X .* w
_scale_rows(X::SparseMatrixCSC{Float64,Int}, w::Vector{Float64}) =
    SparseMatrixCSC{Float64,Int}(Diagonal(w) * X)

"""
Weighted embedding of duplicated compositions (Eq. 5). Returns the weighted matrix and
the indices of the (non-empty) rows that were kept.
"""
function _weight_composition(X::M, scheme::Symbol) where {M<:DenseOrSparse}
    keep = findall(_nonzero_rows(X))
    X = X[keep, :]
    isempty(keep) && return (X, keep)
    g = _row_groups(X)
    n_k = Float64.(g.counts[g.group])
    N = length(keep)
    m = length(g.counts)
    w = scheme === :paper   ? (m .* n_k) ./ N :
        scheme === :bounded ? n_k ./ maximum(g.counts) :
                              (N .* n_k) ./ sum(abs2, g.counts)
    return (_scale_rows(X, w)::M, keep)
end
