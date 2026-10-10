"""
    compressed_angle(S; m = 25_000 * size(S, 1), rng = Random.default_rng()) -> Float64

Hypervolume rescaling factor `Ω(S)^(1/d)` induced by species similarity or functional
complementarity (Song et al. 2025, Eqs. 8–9; R: `get_compressed_angle`). With a
similarity matrix `S` the axis of species `j` points along the `j`-th column of `S`, and
the hypervolume is rescaled by the normalised solid angle of the cone spanned by the
columns of `S`:

    Ω(S) = 2^d |det S| (2π)^(-d/2) ∫_{x ≥ 0} exp(-x'S'Sx/2) dx = 2^d P(Z ≥ 0),
    Z ~ N(0, (S'S)⁻¹),

so that `β_vol = d (Ω(S) vol(P))^(1/d) = Ω(S)^(1/d) ⋅ d vol(P)^(1/d)`. Pass the result
as `betavolume(Z; compressed_angle = compressed_angle(S))`.

- `S = I` (totally dissimilar species): 1.
- Positive off-diagonal similarities shrink the hypervolume (< 1); a singular `S`
  (reciprocal 1-norm condition number below `sqrt(eps)`, e.g. a matrix of ones) gives 0.
- Negative off-diagonal entries (complementarity) expand it (> 1).

With `ρ` the correlation matrix of `(S'S)⁻¹`, the orthant probability is exact for
`d ≤ 3`: `p = 1/4 + asin(ρ₁₂)/(2π)` (d = 2) and
`p = 1/8 + (asin ρ₁₂ + asin ρ₁₃ + asin ρ₂₃)/(4π)` (d = 3). For `d > 3` it is estimated
with the randomised quasi-Monte Carlo algorithm of Genz (MvNormalCDF.jl) using `m`
integration points and the random number generator `rng` (pass a seeded RNG for
reproducible results). `S` need not be symmetric.
"""
function compressed_angle(S::AbstractMatrix{<:Real}; m::Integer = 25_000 * size(S, 1),
                          rng::AbstractRNG = Random.default_rng())
    (size(S, 1) == size(S, 2) && size(S, 1) > 0) ||
        throw(ArgumentError("`S` must be a non-empty square numeric matrix."))
    A = Matrix{Float64}(S)
    _check_finite(A, "S")
    d = size(A, 1)
    _rcond1(A) < sqrt(eps()) && return 0.0
    d == 1 && return 1.0
    p = _solid_angle_fraction(A, m, rng)
    return (2.0^d * max(p, 0.0))^(1 / d)
end

compressed_angle(S; kwargs...) =
    throw(ArgumentError("`S` must be a non-empty square numeric matrix; got a $(typeof(S))."))

"""
    get_compressed_angle(S; kwargs...)

Alias of [`compressed_angle`](@ref), with the name of the R function.
"""
const get_compressed_angle = compressed_angle

"Reciprocal condition number in the 1-norm (LAPACK dgecon, as R's `rcond`)."
function _rcond1(A::Matrix{Float64})
    F = lu(A; check = false)
    issuccess(F) || return 0.0
    anorm = opnorm(A, 1)
    anorm == 0 && return 0.0
    r = try
        LinearAlgebra.LAPACK.gecon!('1', copy(F.factors), anorm)
    catch
        1 / (anorm * opnorm(inv(F), 1))
    end
    return isfinite(r) ? r : 0.0
end

"Fraction of R^d covered by the cone spanned by the columns of S: P(Z ≥ 0)."
function _solid_angle_fraction(S::Matrix{Float64}, m::Integer, rng::AbstractRNG)
    d = size(S, 1)
    Σ = inv(S' * S)
    Σ = (Σ .+ Σ') ./ 2
    is = sqrt.(1 ./ diag(Σ))
    ρ = similar(Σ)
    for j in 1:d, i in 1:d
        ρ[i, j] = i == j ? 1.0 : clamp((is[i] * Σ[i, j]) * is[j], -1.0, 1.0)
    end
    if d == 2
        return 1 / 4 + asin(ρ[1, 2]) / (2π)
    elseif d == 3
        return 1 / 8 + (asin(ρ[1, 2]) + asin(ρ[1, 3]) + asin(ρ[2, 3])) / (4π)
    end
    p, _ = MvNormalCDF.mvnormcdf(zeros(d), ρ, zeros(d), fill(Inf, d); m = m, rng = rng)
    return p
end
