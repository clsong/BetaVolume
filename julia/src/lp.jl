# A small dense simplex method, used to find a strictly interior point (the Chebyshev
# centre) of the intersection of two convex hulls.

"""
    _simplex_max(c, M, b; maxiter, tol) -> (status, z, objective)

Solve `max c'z  s.t.  M z ≤ b, z ≥ 0` for `b ≥ 0` (so `z = 0` is feasible and no phase 1
is needed) with the condensed (Tucker) tableau and Bland's anti-cycling rule. `status`
is `:optimal`, `:unbounded` or `:iteration_limit`; in the last case `z` is the current
(feasible) basic solution.
"""
function _simplex_max(c::AbstractVector{Float64}, M::AbstractMatrix{Float64},
                      b::AbstractVector{Float64};
                      maxiter::Integer = 50 * (size(M, 1) + size(M, 2)) + 1000,
                      tol::Float64 = 1e-12)
    m, n = size(M)
    length(c) == n && length(b) == m || throw(DimensionMismatch("inconsistent LP sizes"))
    all(>=(0), b) || throw(ArgumentError("_simplex_max requires b ≥ 0"))
    T = Matrix{Float64}(undef, m + 1, n + 1)
    T[1:m, 1:n] .= M
    T[1:m, n + 1] .= b
    T[m + 1, 1:n] .= .-c       # objective row stores -c (f = rhs - Σ T x_N)
    T[m + 1, n + 1] = 0.0
    col_lab = collect(1:n)            # labels of nonbasic variables (columns)
    row_lab = collect((n + 1):(n + m))  # labels of basic variables (rows); slacks > n
    status = :iteration_limit
    @inbounds for _ in 1:maxiter
        # entering variable: smallest label with negative entry in the objective row
        e = 0
        best = typemax(Int)
        for j in 1:n
            if T[m + 1, j] < -tol && col_lab[j] < best
                e = j
                best = col_lab[j]
            end
        end
        if e == 0
            status = :optimal
            break
        end
        # ratio test, ties broken by the smallest label
        p = 0
        minratio = Inf
        plab = typemax(Int)
        for i in 1:m
            a = T[i, e]
            a > tol || continue
            ratio = T[i, n + 1] / a
            tie = abs(ratio - minratio) <= 1e-13 * max(1.0, abs(minratio))
            if (ratio < minratio && !tie) || (tie && row_lab[i] < plab)
                p = i
                minratio = min(ratio, minratio)
                plab = row_lab[i]
            end
        end
        if p == 0
            status = :unbounded
            break
        end
        # pivot on (p, e)
        piv = T[p, e]
        for j in 1:(n + 1)
            j == e && continue
            T[p, j] /= piv
        end
        for j in 1:(n + 1)
            j == e && continue
            tpj = T[p, j]
            tpj == 0 && continue
            for i in 1:(m + 1)
                i == p && continue
                T[i, j] -= T[i, e] * tpj
            end
        end
        for i in 1:(m + 1)
            i == p && continue
            T[i, e] = -T[i, e] / piv
        end
        T[p, e] = 1 / piv
        for i in 1:m   # remove round-off infeasibility
            if T[i, n + 1] < 0
                T[i, n + 1] = 0.0
            end
        end
        col_lab[e], row_lab[p] = row_lab[p], col_lab[e]
    end
    z = zeros(n)
    for i in 1:m
        row_lab[i] <= n && (z[row_lab[i]] = T[i, n + 1])
    end
    return (status, z, dot(c, z))
end

"""
    _chebyshev_center(A, b) -> (x, r, status)

Centre `x` and radius `r` of the largest ball inside `{x : A x ≤ b}`, i.e.
`max r  s.t.  a_i'x + ‖a_i‖ r ≤ b_i`, for `b ≥ 0` (the origin is feasible; tiny negative
entries of `b` are clamped to 0). The free vector `x` is split as `x = u - v`.
"""
function _chebyshev_center(A::AbstractMatrix{Float64}, b::AbstractVector{Float64}; kwargs...)
    f, d = size(A)
    nrm = [norm(view(A, i, :)) for i in 1:f]
    M = hcat(A, -A, nrm)
    c = zeros(2d + 1)
    c[end] = 1.0
    status, z, _ = _simplex_max(c, M, max.(b, 0.0); kwargs...)
    x = z[1:d] .- z[(d + 1):(2d)]
    r = status === :unbounded ? Inf : z[end]
    return (x, r, status)
end
