"""
    contributions(Z; weights = false, remove_unique = false, method = :auto,
                  dim_threshold = 10, preprocess = true, weight_scheme = :legacy)
        -> (values, units)

Contribution of each embedded point (a community or a species, whichever is the more
numerous dimension of `Z`) to geometric beta diversity, as the normalised loss of
beta diversity when the point is removed (Song et al. 2025, Eq. 7 and Fig. 6;
R: `disentangle_individual_contribution`):

    C_i = (vol(P₀)^(1/d) − vol(P₋ᵢ)^(1/d)) / Σ_j (vol(P₀)^(1/d) − vol(P₋ⱼ)^(1/d)).

# Units of removal
Without weighting each row of the embedding is removed in turn; if a composition occurs
several times, removing one copy leaves `P` unchanged, so every copy contributes 0. With
`weights = true` (and `preprocess = true`), identical rows of the weighted embedding are
merged into one unit and each *unique composition* is removed in turn, as in Fig. 6. The
weights are computed once on the full metacommunity. `:paper` and `:legacy` weights give
identical contributions.

# Methods
- exact: only hull vertices occurring once are recomputed (other points contribute
  exactly 0); losses are clamped at 0, so contributions are non-negative.
- hyper-normal: all leave-one-unit-out determinants come from one Cholesky factorisation
  via the matrix determinant lemma. With `n` points (origin included), mean `μ`,
  covariance `C` and a unit of `n_k` identical points `x_k`, `δ = x_k − μ`:

      log det C₋ₖ = log det C + d log(n−1) + log(1 − n n_k/((n−n_k)(n−1)) δ'C⁻¹δ) − d log(n−n_k−1).

  Removing a point near the centroid can increase the determinant, so small negative
  contributions are possible.

If beta diversity is zero (degenerate embedding) all values are `NaN` and a warning is
issued. If no single unit matters (every hull vertex occurs more than once, without
weights), all values are 0 and a warning suggests `weights = true`.

# Returns
A `NamedTuple` with `values::Vector{Float64}` (summing to one) and
`units::Vector{Vector{Int}}`, where `units[k]` lists the `ids` (see [`preprocess`](@ref):
rows of `Z`, or columns when `Z` was transposed) forming unit `k`.

# Example
```julia
fig4 = [1 0; 0 1; 0 1; 1 1; 1 1; 1 1]
c = contributions(fig4; weights = true)
c.values  # ≈ [0.1617, 0.3725, 0.4658]
c.units   # [[1], [2, 3], [4, 5, 6]]
```
"""
function contributions(Z; weights::Bool = false, remove_unique::Bool = false,
                       method = :auto, dim_threshold::Real = 10, preprocess::Bool = true,
                       weight_scheme = :legacy)
    scheme = _resolve_scheme(weight_scheme)
    if preprocess
        E = _preprocess(Z, weights, remove_unique, scheme)
        X = E.points
        ids::Vector{Int} = E.ids
    else
        X = _as_matrix(Z, "Z")
        ids = collect(1:size(X, 1))
    end
    n, d = size(X)
    m = _resolve_method(method, d, dim_threshold)
    n == 0 && return (values = Float64[], units = Vector{Int}[])

    grouped = weights && preprocess
    if grouped
        g = _row_groups(X)
        group, first = g.group, g.first
        units = Vector{Int}[Int[] for _ in first]
        for i in 1:n
            push!(units[group[i]], ids[i])
        end
    else
        group = collect(1:n)
        first = group
        units = Vector{Int}[[ids[i]] for i in 1:n]
    end
    nunits = length(units)

    loss, scale, degenerate = if m === :exact
        Xd = Matrix(X)
        _exact_loo_loss(grouped ? Xd[first, :] : Xd)
    else
        _hn_loo_loss(X, group, first)
    end

    if degenerate
        @warn "The metacommunity has zero beta diversity; contributions are undefined."
        return (values = fill(NaN, nunits), units = units)
    end
    total = sum(loss)
    if !isfinite(total) || total <= 1e-10 * scale
        @warn "Removing any single point leaves beta diversity unchanged (every vertex of " *
              "the hull occurs more than once), so all contributions are 0. Use " *
              "`weights = true` to treat duplicated compositions as single units."
        return (values = zeros(nunits), units = units)
    end
    return (values = loss ./ total, units = units)
end

"""
    disentangle_individual_contribution(Z; kwargs...)

Alias of [`contributions`](@ref), with the name of the R function.
"""
const disentangle_individual_contribution = contributions

"""
Leave-one-out loss of vol^(1/d) for the exact hull. Returns `(loss, scale, degenerate)`.
Only hull vertices whose composition occurs once can change the hull when removed.
"""
function _exact_loo_loss(X::Matrix{Float64})
    n, d = size(X)
    loss = zeros(n)
    g = _row_groups(X)
    keepg = [k for k in eachindex(g.first) if any(!iszero, view(X, g.first[k], :))]
    owner = g.first[keepg]          # row represented by each unique non-origin point
    mult = g.counts[keepg]
    nu = length(owner)
    U = zeros(d, nu + 1)            # unique points, origin last
    for (k, i) in enumerate(owner)
        U[:, k] .= view(X, i, :)
    end
    h = _hull(U; vertices = true)
    h.degenerate && return (loss, NaN, true)
    v0 = h.volume^(1 / d)
    cands = [k for k in h.vertices if k <= nu && mult[k] == 1]
    vk = zeros(length(cands))
    loo(j) = begin
        k = cands[j]
        Uk = U[:, [i for i in 1:(nu + 1) if i != k]]
        vk[j] = _hull_volume_cols(Uk)^(1 / d)
    end
    if Threads.nthreads() > 1 && length(cands) > 1
        Threads.@threads for j in eachindex(cands)
            loo(j)
        end
    else
        foreach(loo, eachindex(cands))
    end
    for (j, k) in enumerate(cands)
        loss[owner[k]] = max(v0 - vk[j], 0.0)
    end
    return (loss, v0, false)
end

"""
Leave-one-group-out loss of det(C)^(1/d) for the hyper-normal approximation (the factor
4d cancels in the normalised contributions). Returns `(loss, scale, degenerate)`.
"""
function _hn_loo_loss(X::DenseOrSparse, group::Vector{Int}, first::Vector{Int})
    n0, d = size(X)
    n = n0 + 1
    nunits = length(first)
    n_k = zeros(Int, nunits)
    for gi in group
        n_k[gi] += 1
    end
    C = _embedding_cov(X)
    ld, F = _logdet_psd(C)
    (F === nothing || !isfinite(ld)) && return (zeros(nunits), NaN, true)
    μ = vec(sum(X; dims = 1)) ./ n
    Δ = Matrix(X[first, :])' .- μ                 # d × nunits
    Zs = LowerTriangular(F.U') \ Δ                 # C = U'U
    maha = vec(sum(abs2, Zs; dims = 1))
    beta0 = exp(ld / d)
    logdet_scatter = ld + d * log(n - 1)
    loss = Vector{Float64}(undef, nunits)
    for k in 1:nunits
        nk = Float64(n_k[k])
        remaining = n - nk
        factor = 1 - n * nk / ((n - nk) * (n - 1)) * maha[k]
        bl = (remaining >= 2 && factor > sqrt(eps())) ?
             exp((logdet_scatter + log(factor) - d * log(remaining - 1)) / d) : 0.0
        loss[k] = beta0 - bl
    end
    return (loss, beta0, false)
end
