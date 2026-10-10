"""
    temporal_overlap(Z1, Z2; oriented = false, preprocess = true, weights = false,
                     weight_scheme = :legacy, remove_unique = false) -> Float64

Temporal change of beta diversity as the (oriented) overlap of the hypervolumes at two
times (Song et al. 2025, Eq. 6 and Fig. 5; R: `beta_temporal_overlap`):

    β_{t→t+1} = (vol(P_{t+1} ∩ P_t) / vol(P_t))^(1/d).

- **Unoriented** (`oriented = false`): ratio of the volume of the intersection of the two
  convex hulls (both including the origin) to the volume of the hull at time 1, in
  [0, 1]: 1 if `P₂` covers `P₁`, 0 if the hulls only share a lower-dimensional face.
- **Oriented** (`oriented = true`): requires each embedding to be a simplex, i.e. exactly
  `d` embedded points, in the same order at both times. With signed volumes
  `s_i = det(X_i)/d!`: if the orientations agree (or `s₂ = 0`) the overlap is
  `vol(P₁ ∩ P₂)`, otherwise `|s₁| + |s₂|`; the result is `(overlap/|s₁|)^(1/d)`. Values
  above 1 indicate synchronous changes in opposite directions.

`Z1` and `Z2` must have the same dimensions (the same species and communities in the same
order). `weights`/`weight_scheme` are passed to [`preprocess`](@ref);
`remove_unique = true` is rejected because it can misalign the axes of the two
embeddings. Returns `NaN` (with a warning) when the metacommunity at time 1 has zero
volume.

The intersection volume is computed from the facet inequalities of both hulls: a strictly
interior point (Chebyshev centre) is found with a small simplex method, the vertices of
the intersection follow from the dual convex hull (halfspace intersection) and their
convex hull gives the volume.

# Example (Fig. 5)
```julia
t0 = [0 1; 1 1]; t1 = [0 1; 1 0]; t2 = [1 0; 0 1]
temporal_overlap(t0, t1)                   # √(1/2)
temporal_overlap(t1, t2)                   # 1
temporal_overlap(t1, t2; oriented = true)  # √2
```
"""
function temporal_overlap(Z1, Z2; oriented::Bool = false, preprocess::Bool = true,
                          weights::Bool = false, weight_scheme = :legacy,
                          remove_unique::Bool = false)
    scheme = _resolve_scheme(weight_scheme)
    m1 = _as_matrix(Z1, "Z1")
    m2 = _as_matrix(Z2, "Z2")
    size(m1) == size(m2) || throw(ArgumentError(
        "`Z1` and `Z2` must have the same dimensions (the same species and communities " *
        "in the same order)."))
    if preprocess
        remove_unique && throw(ArgumentError(
            "`remove_unique = true` is not supported by temporal_overlap() because it can " *
            "misalign the axes of the two embeddings."))
        m1 = _preprocess_matrix(m1, weights, false, scheme).points
        m2 = _preprocess_matrix(m2, weights, false, scheme).points
    end
    X1 = Matrix(m1)
    X2 = Matrix(m2)
    size(X1, 2) == size(X2, 2) ||
        throw(ArgumentError("The two embeddings have different dimensions."))
    d = size(X1, 2)
    P1 = _hull_points(X1)
    P2 = _hull_points(X2)

    if !oriented
        v1 = _hull_volume_cols(P1)
        if v1 <= 0
            @warn "The metacommunity at time 1 has zero volume; the overlap is undefined."
            return NaN
        end
        v12 = _intersection_volume(P1, P2)
        return min(v12 / v1, 1.0)^(1 / d)
    end

    (size(X1, 1) == d && size(X2, 1) == d) || throw(ArgumentError(
        "The oriented overlap requires each embedding to be a simplex: exactly $d " *
        "embedded points (non-empty communities) in $d dimensions."))
    s1 = det(X1) / _factorial(d)
    s2 = det(X2) / _factorial(d)
    if s1 == 0
        @warn "The metacommunity at time 1 has zero volume; the overlap is undefined."
        return NaN
    end
    overlap = (s2 == 0 || sign(s1) == sign(s2)) ? _intersection_volume(P1, P2) :
              abs(s1) + abs(s2)
    return (overlap / abs(s1))^(1 / d)
end

"""
    beta_temporal_overlap(Z1, Z2; kwargs...)

Alias of [`temporal_overlap`](@ref), with the name of the R function.
"""
const beta_temporal_overlap = temporal_overlap

"""
    _intersection_volume(P1, P2)

Volume of the intersection of the convex hulls of the columns of `P1` and `P2` (d × m,
distinct points, both containing the origin).
"""
function _intersection_volume(P1::Matrix{Float64}, P2::Matrix{Float64})
    d = size(P1, 1)
    if d == 1
        lo = max(minimum(P1), minimum(P2))
        hi = min(maximum(P1), maximum(P2))
        return max(hi - lo, 0.0)
    end
    (_is_full_dimensional(P1) && _is_full_dimensional(P2)) || return 0.0
    h1 = _hull(P1; equations = true)
    h2 = _hull(P2; equations = true)
    (h1.degenerate || h2.degenerate) && return 0.0
    E = hcat(h1.equations, h2.equations)
    E = E[:, _col_groups(E).first]                 # drop duplicated (coplanar) facets
    scale = max(maximum(abs, P1), maximum(abs, P2))
    A = permutedims(E[1:d, :])
    b = -E[d + 1, :] ./ scale                      # origin is a vertex of both: b ≥ 0
    x, r, _ = _chebyshev_center(A, b)
    (isfinite(r) && r > 1e-10) || return 0.0       # lower-dimensional intersection
    V = _halfspace_vertices(A, -E[d + 1, :], x .* scale)
    V === nothing && return 0.0
    return _hull_volume_cols(_hull_points(permutedims(V)))
end

"""
    _halfspace_vertices(A, b, x0) -> d × k matrix or `nothing`

Vertices of the bounded polytope `{x : A x ≤ b}` given a strictly interior point `x0`,
by duality (as Qhull's halfspace mode): with `β_i = b_i − a_i'x0 > 0` the dual points are
`a_i / β_i`, and each facet `n'y + c ≤ 0` of their convex hull gives the vertex
`x0 − n / c`.
"""
function _halfspace_vertices(A::Matrix{Float64}, b::Vector{Float64}, x0::Vector{Float64})
    f, d = size(A)
    β = b .- A * x0
    all(>(0), β) || return nothing
    D = permutedims(A ./ β)                         # d × f dual points
    D = D[:, _col_groups(D).first]
    D = hcat(D, zeros(d))                           # origin is interior: harmless for the hull
    h = _hull(D; equations = true)
    h.degenerate && return nothing
    Eq = h.equations
    V = Eq[1:d, :] ./ (.-Eq[d + 1:d + 1, :]) .+ x0
    return V[:, _col_groups(V).first]
end
