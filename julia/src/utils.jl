# Internal helpers: input validation, grouping of identical rows, rank, argument
# resolution. Nothing in this file is exported.

const DenseOrSparse = Union{Matrix{Float64},SparseMatrixCSC{Float64,Int}}

function _check_finite(values, arg::AbstractString)
    any(isnan, values) && throw(ArgumentError("`$arg` contains missing values (NaN)."))
    all(isfinite, values) || throw(ArgumentError("`$arg` contains infinite values."))
    return nothing
end

function _check_nonempty(Z, arg::AbstractString)
    (size(Z, 1) == 0 || size(Z, 2) == 0) &&
        throw(ArgumentError("`$arg` must have at least one row and one column."))
    return nothing
end

const _SparseLike = Union{SparseArrays.AbstractSparseMatrix{<:Real},
                          Adjoint{<:Real,<:SparseArrays.AbstractSparseMatrix},
                          Transpose{<:Real,<:SparseArrays.AbstractSparseMatrix}}

"""
    _as_matrix(Z, arg) -> Matrix{Float64} or SparseMatrixCSC{Float64,Int}

Validate user input and convert it to a fresh (never aliased) Float64 matrix. Sparse
inputs stay sparse (with explicitly stored zeros dropped).
"""
function _as_matrix(Z::_SparseLike, arg::AbstractString)
    _check_nonempty(Z, arg)
    X = copy(SparseMatrixCSC{Float64,Int}(sparse(Z)))
    _check_finite(nonzeros(X), arg)
    dropzeros!(X)
    return X
end

function _as_matrix(Z::AbstractMatrix{<:Real}, arg::AbstractString)
    _check_nonempty(Z, arg)
    X = Matrix{Float64}(Z)  # always a copy
    _check_finite(X, arg)
    return X
end

_as_matrix(Z::AbstractVector{<:Real}, arg::AbstractString) =
    _as_matrix(reshape(collect(Z), 1, :), arg)

function _as_matrix(Z::AbstractMatrix, arg::AbstractString)
    any(ismissing, Z) && throw(ArgumentError("`$arg` contains missing values."))
    all(x -> x isa Real, Z) ||
        throw(ArgumentError("`$arg` must contain numeric (real) values."))
    return _as_matrix(Float64.(Z), arg)
end

_as_matrix(Z, arg::AbstractString) = throw(ArgumentError(
    "`$arg` must be a numeric matrix (AbstractMatrix{<:Real}); got a $(typeof(Z))."))

# --- grouping of identical rows ---------------------------------------------------

"""
    _row_groups(X) -> (group, counts, first)

Group identical rows (exact equality; -0.0 == 0.0). `group[i]` is the id of the unique
row of row `i` (ids in order of first appearance), `counts` the multiplicities and
`first` the index of the first row of each group.
"""
function _row_groups(X::Matrix{Float64})
    n, d = size(X)
    h = zeros(UInt, n)
    @inbounds for j in 1:d, i in 1:n
        h[i] = hash(X[i, j] + 0.0, h[i])
    end
    rows_equal(i, k) = all(j -> @inbounds(X[i, j] == X[k, j]), 1:d)
    return _group_by_hash(h, rows_equal)
end

_row_groups(X::SparseMatrixCSC{Float64,Int}) = _col_groups(permutedims(X))

function _col_groups(X::SparseMatrixCSC{Float64,Int})
    n = size(X, 2)
    rv = rowvals(X)
    nz = nonzeros(X)
    h = zeros(UInt, n)
    @inbounds for j in 1:n
        hj = UInt(0x5bd1e995)
        for p in nzrange(X, j)
            nz[p] == 0 && continue
            hj = hash(rv[p], hash(nz[p] + 0.0, hj))
        end
        h[j] = hj
    end
    function cols_equal(i, k)
        ri = [(rv[p], nz[p]) for p in nzrange(X, i) if nz[p] != 0]
        rk = [(rv[p], nz[p]) for p in nzrange(X, k) if nz[p] != 0]
        length(ri) == length(rk) || return false
        return all(t -> t[1][1] == t[2][1] && t[1][2] == t[2][2], zip(ri, rk))
    end
    return _group_by_hash(h, cols_equal)
end

_col_groups(X::Matrix{Float64}) = _row_groups(permutedims(X))

function _group_by_hash(h::Vector{UInt}, equal::F) where {F}
    n = length(h)
    group = Vector{Int}(undef, n)
    first = Int[]
    counts = Int[]
    seen = Dict{UInt,Vector{Int}}()
    for i in 1:n
        gs = get!(() -> Int[], seen, h[i])
        found = 0
        for g in gs
            if equal(i, first[g])
                found = g
                break
            end
        end
        if found == 0
            push!(first, i)
            push!(counts, 1)
            found = length(first)
            push!(gs, found)
        else
            counts[found] += 1
        end
        group[i] = found
    end
    return (group = group, counts = counts, first = first)
end

# --- rank ------------------------------------------------------------------------

"""
    _qr_rank!(W, tol)

Rank of `W` (n × p, overwritten) with the criterion of LINPACK `dqrdc2`, i.e. of R's
`qr(x, tol)$rank`: columns are processed in order and a column is declared negligible
when its norm, after projecting out the previously accepted columns, is below `tol`
times its original norm. Householder reflections are used.
"""
function _qr_rank!(W::Matrix{Float64}, tol::Float64)
    n, p = size(W)
    rank = 0
    l = 1
    @inbounds for j in 1:p
        l > n && break
        orig = norm(view(W, :, j))   # invariant under the reflections applied so far
        orig == 0 && (orig = 1.0)
        x = view(W, l:n, j)
        r = norm(x)
        r < tol * orig && continue
        rank += 1
        if l < n
            α = x[1] > 0 ? -r : r
            v1 = x[1] - α
            vtv = v1^2
            for i in 2:length(x)
                vtv += x[i]^2
            end
            for k in (j + 1):p
                y = view(W, l:n, k)
                s = v1 * y[1]
                for i in 2:length(x)
                    s += x[i] * y[i]
                end
                c = 2s / vtv
                y[1] -= c * v1
                for i in 2:length(x)
                    y[i] -= c * x[i]
                end
            end
            x[1] = α
            for i in 2:length(x)
                x[i] = 0.0
            end
        end
        l += 1
    end
    return min(rank, n)
end

"""
    _is_full_dimensional(P; tol = 1e-10)

`P` is `d × m` with the points as columns, the origin among them. The convex hull is
full dimensional iff the points span R^d linearly (R: `is_full_dimensional`).
"""
function _is_full_dimensional(P::AbstractMatrix{Float64}; tol::Float64 = 1e-10)
    d, m = size(P)
    m < d + 1 && return false
    scale = maximum(abs, P)
    scale == 0 && return false
    W = permutedims(P) ./ scale
    return _qr_rank!(W, tol) >= d
end

# --- argument resolution -----------------------------------------------------------

function _resolve_method(method, d::Integer, dim_threshold::Real)
    m = method isa AbstractString ? Symbol(method) : method
    m === nothing && (m = :auto)
    if m === :auto
        return d > dim_threshold ? :hyper_normal : :exact
    elseif m === :exact || m === :deterministic
        return :exact
    elseif m === :hyper_normal
        return :hyper_normal
    end
    throw(ArgumentError("`method` must be one of :auto, :exact, :deterministic or " *
                        ":hyper_normal; got $(repr(method))."))
end

function _resolve_scheme(weight_scheme)
    s = weight_scheme isa AbstractString ? Symbol(weight_scheme) : weight_scheme
    (s === :paper || s === :legacy) && return s
    throw(ArgumentError("`weight_scheme` must be :paper or :legacy; got $(repr(weight_scheme))."))
end

_factorial(d::Integer) = prod(Float64, 1:d; init = 1.0)
