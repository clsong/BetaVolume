# Thin, checked interface to Qhull.
#
# `DirectQhull.ConvexHull` (v0.2) ignores Qhull's exit code, prints errors to stderr
# and, for every call, converts the whole facet list to Julia arrays, which costs more
# than the hull itself in the leave-one-out loops. We therefore call `qh_new_qhull` of
# the reentrant Qhull library shipped with DirectQhull (Qhull_jll) directly:
#
# * the exit code is checked (2 = qh_ERRsingular, i.e. flat input);
# * Qhull's diagnostic output goes to the null device instead of stderr;
# * the volume is read with the library getters (`qh_getarea`, `qh_get_totvol`), and
#   vertices / facet equations are extracted only when needed, with DirectQhull's own
#   (internal) helpers.
#
# The options are those of the R package: "Qt" (triangulated output, no "Tv").
# Qhull itself adds "Qx" (exact pre-merges) in 5-d and higher.

const _DQ = DirectQhull
const _QH_ERRFILE = Ref{Ptr{Cvoid}}(C_NULL)
const _QH_OPTIONS = "Qt"

function _init_qhull_errfile()
    try
        io = open(Sys.iswindows() ? "NUL" : "/dev/null", "w")
        f = Libc.FILE(io)        # duplicates the descriptor
        close(io)
        _QH_ERRFILE[] = f.ptr
    catch
        _QH_ERRFILE[] = C_NULL   # fall back to stderr
    end
    return nothing
end

"""
    _qhull(f, P, opts)

Run Qhull on the points in the columns of `P` (d × m) and call `f(qh, code)` with the
Qhull pointer and the exit code before freeing all Qhull memory.
"""
function _qhull(f::F, P::Matrix{Float64}, opts::String) where {F}
    qh = _DQ.qhullptr()
    try
        code = ccall((:qh_new_qhull, _DQ.qh_lib), Cint,
                     (Ptr{_DQ.qhT}, Cint, Cint, Ptr{Float64}, Cuint, Cstring, Ptr{Cvoid},
                      Ptr{Cvoid}),
                     qh, size(P, 1), size(P, 2), P, Cuint(0), "qhull " * opts, C_NULL,
                     _QH_ERRFILE[])
        return f(qh, Int(code))
    finally
        _DQ.qhullfree(qh)
    end
end

function _qh_totvol(qh)
    fl = ccall((:qh_get_facet_list, _DQ.qh_lib), Ptr{Cvoid}, (Ptr{_DQ.qhT},), qh)
    ccall((:qh_getarea, _DQ.qh_lib), Cvoid, (Ptr{_DQ.qhT}, Ptr{Cvoid}), qh, fl)
    return ccall((:qh_get_totvol, _DQ.qh_lib), Cdouble, (Ptr{_DQ.qhT},), qh)
end

const _EMPTY_EQ = Matrix{Float64}(undef, 0, 0)

function _qhull_run(P::Matrix{Float64}, opts::String, want_vertices::Bool,
                    want_equations::Bool)
    d = size(P, 1)
    return _qhull(P, opts) do qh, code
        code == 0 || return (code, 0.0, Int[], _EMPTY_EQ)
        vol = _qh_totvol(qh)
        verts = want_vertices ?
                sort!(unique!(Int.(vec(_DQ.qh_get_convex_hull_pnts(qh, Val(d)))))) : Int[]
        eqs = want_equations ?
              Matrix{Float64}(_DQ.qh_get_simplex_facet_arrays(qh, Val(d))[3]) : _EMPTY_EQ
        return (0, vol, verts, eqs)
    end
end

"""
    _hull(P; vertices = false, equations = false)

Convex hull of the distinct points in the columns of `P` (d × m), which must contain
the origin. Returns `(volume, vertices, equations, degenerate)`; a hull that is not full
dimensional has volume 0 and `degenerate = true` (Qhull is never called on such input).
`vertices` are column indices of `P`; `equations` is `(d+1) × nfacets` with unit outward
normals and offsets (`n'x + offset ≤ 0` inside).
"""
function _hull(P::Matrix{Float64}; vertices::Bool = false, equations::Bool = false)
    d, m = size(P)
    if d == 1
        lo, hi = extrema(view(P, 1, :))
        vol = max(hi, 0.0) - min(lo, 0.0)
        vol <= 0 && return (volume = 0.0, vertices = Int[], equations = _EMPTY_EQ,
                            degenerate = true)
        return (volume = vol, vertices = unique([argmin(view(P, 1, :)), argmax(view(P, 1, :))]),
                equations = _EMPTY_EQ, degenerate = false)
    end
    if !_is_full_dimensional(P)
        return (volume = 0.0, vertices = Int[], equations = _EMPTY_EQ, degenerate = true)
    end
    code, vol, verts, eqs = _qhull_run(P, _QH_OPTIONS, vertices, equations)
    if code == 1 || code == 2
        # input / singular error: precision-level flat input, zero volume (as in R)
        return (volume = 0.0, vertices = Int[], equations = _EMPTY_EQ, degenerate = true)
    elseif code != 0
        # precision or topology error: retry once with joggled input (as in R)
        code2, vol, verts, eqs = _qhull_run(P, _QH_OPTIONS * " QJ", vertices, equations)
        code2 == 0 || error("Qhull failed to compute the convex hull (exit code $code).")
    end
    return (volume = vol, vertices = verts, equations = eqs, degenerate = false)
end

_hull_volume_cols(P::Matrix{Float64}) = _hull(P).volume

"""
    _hull_points(X) -> P

Distinct rows of `X` (n × d) plus the origin, as the columns of a `d × m` matrix.
"""
function _hull_points(X::AbstractMatrix{Float64})
    n, d = size(X)
    Y = Matrix{Float64}(undef, n + 1, d)
    Y[1:n, :] .= X
    Y[n + 1, :] .= 0.0
    g = _row_groups(Y)
    return permutedims(Y[g.first, :])
end

"""
    _hull_volume(X)

Exact volume of the convex hull of the rows of `X` together with the origin.
"""
_hull_volume(X::AbstractMatrix{Float64}) = _hull_volume_cols(_hull_points(X))
