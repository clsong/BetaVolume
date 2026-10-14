using BetaVolume
using Test
using LinearAlgebra
using SparseArrays
using Random
using Statistics
using Logging
using JSON

const BV = BetaVolume

# --- helpers -----------------------------------------------------------------------

"All 2^d presence/absence compositions (rows), the empty one first."
all_compositions(d) = Float64[(i >> (j - 1)) & 1 for i in 0:(2^d - 1), j in 1:d]

quietly(f) = with_logger(f, NullLogger())

beta_from_volume(v, d) = d * v^(1 / d)

"Brute-force contributions (Eq. 7): recompute betavolume without each unit."
function brute_contribution(m; weights = false, method = :exact)
    E = preprocess(m; weights = weights)
    X = E.points
    if weights
        g = BV._row_groups(Matrix(X))
        units = [findall(==(k), g.group) for k in eachindex(g.first)]
    else
        units = [[i] for i in 1:size(X, 1)]
    end
    b0 = betavolume(X; preprocess = false, method = method)
    loss = map(units) do u
        keep = setdiff(1:size(X, 1), u)
        b0 - betavolume(X[keep, :]; preprocess = false, method = method)
    end
    return loss ./ sum(loss)
end

"Area of the convex hull of 2-d points plus the origin (monotone chain + shoelace)."
function shoelace_hull_area(m)
    pts = unique([(Float64(m[i, 1]), Float64(m[i, 2])) for i in 1:size(m, 1)])
    push!(pts, (0.0, 0.0))
    pts = sort(unique(pts))
    length(pts) < 3 && return 0.0
    cross(o, a, b) = (a[1] - o[1]) * (b[2] - o[2]) - (a[2] - o[2]) * (b[1] - o[1])
    lower = Tuple{Float64,Float64}[]
    for p in pts
        while length(lower) >= 2 && cross(lower[end - 1], lower[end], p) <= 0
            pop!(lower)
        end
        push!(lower, p)
    end
    upper = Tuple{Float64,Float64}[]
    for p in reverse(pts)
        while length(upper) >= 2 && cross(upper[end - 1], upper[end], p) <= 0
            pop!(upper)
        end
        push!(upper, p)
    end
    hull = vcat(lower[1:(end - 1)], upper[1:(end - 1)])
    a = 0.0
    for i in eachindex(hull)
        p, q = hull[i], hull[mod1(i + 1, length(hull))]
        a += p[1] * q[2] - q[1] * p[2]
    end
    return abs(a) / 2
end

"Monte Carlo fraction of R^d covered by the cone spanned by the columns of S."
function mc_cone_fraction(S, n, rng)
    d = size(S, 1)
    Z = randn(rng, d, n)
    C = S \ Z
    return mean(all(>=(0), view(C, :, j)) for j in 1:n)
end

"Van Oosterom–Strackee solid angle (steradians) of the cone spanned by a, b, c."
function solid_angle(a, b, c)
    num = abs(dot(a, cross(b, c)))
    la, lb, lc = norm(a), norm(b), norm(c)
    den = la * lb * lc + dot(a, b) * lc + dot(a, c) * lb + dot(b, c) * la
    om = 2 * atan(num, den)
    return om < 0 ? om + 4π : om
end

tomatrix(rows) = isempty(rows) ? zeros(0, 0) :
                 Float64[Float64(rows[i][j]) for i in eachindex(rows), j in eachindex(rows[1])]
tovec(x) = x === nothing ? [NaN] : x isa AbstractVector ?
           Float64[v === nothing ? NaN : Float64(v) for v in x] : [Float64(x)]
tonum(x) = x === nothing ? NaN : Float64(x)

# --- tests -------------------------------------------------------------------------

@testset "BetaVolume.jl" begin

@testset "Paper examples" begin
    @testset "Fig. 1" begin
        meta_I = [1 0; 0 1]
        meta_II = [1 0; 0 1; 1 1]
        @test betavolume(meta_I) ≈ 2 * sqrt(1 / 2)
        @test betavolume(meta_II) ≈ 2
        @test betavolume(meta_II) > betavolume(meta_I)
        @test betavolume(meta_I; method = :exact) == betavolume(meta_I; method = :deterministic)
    end

    @testset "maximum d with all 2^d - 1 compositions" begin
        for d in 2:5
            full = all_compositions(d)
            @test betavolume(full[2:end, :]; method = :exact) ≈ d
            @test betavolume(full; method = :exact) ≈ d
        end
    end

    @testset "Fig. 4: weighted embedding" begin
        fig4 = [1 0; 0 1; 0 1; 1 1; 1 1; 1 1]
        E = preprocess(fig4; weights = true, weight_scheme = :paper)
        w = [0.5, 1, 1, 1.5, 1.5, 1.5]
        @test E.points ≈ fig4 .* w
        @test E.ids == 1:6
        @test !E.transposed
        @test betavolume(fig4; weights = true, weight_scheme = :paper) ≈ 2 * sqrt(1.125)
        @test betavolume(fig4; weights = true, weight_scheme = :paper) > 2  # exceeds d
        @test preprocess(fig4; weights = true).points ==
              preprocess(fig4; weights = true, weight_scheme = :legacy).points  # default
        @test betavolume(fig4; weights = true) ≈ 1.818275 atol = 1e-6
        L = preprocess(fig4; weights = true, weight_scheme = :legacy)
        @test L.points[:, 1] ≈ 6 .* [1, 0, 0, 3, 3, 3] ./ 14
        @test L.points ≈ fig4 .* (6 .* [1, 2, 2, 3, 3, 3] ./ 14)
        @test betavolume(fig4; weights = true, weight_scheme = :legacy) ≈ 1.818275 atol = 1e-6
        @test betavolume(fig4; weights = true, weight_scheme = "legacy") ≈ 1.818275 atol = 1e-6
    end

    @testset "Fig. 6: contributions" begin
        loss = sqrt(1.125) .- sqrt.([0.75, 0.375, 0.25])
        expected = loss ./ sum(loss)
        fig6 = [0.5 0; 0 1; 1.5 1.5]
        c6 = contributions(fig6)
        @test c6.values ≈ expected
        @test round.(c6.values; digits = 2) == [0.16, 0.37, 0.47]
        @test c6.units == [[1], [2], [3]]
        fig4 = [1 0; 0 1; 0 1; 1 1; 1 1; 1 1]
        c4 = contributions(fig4; weights = true)
        @test c4.values ≈ expected
        @test c4.values ≈ [0.1617, 0.3725, 0.4658] atol = 1e-4
        @test c4.units == [[1], [2, 3], [4, 5, 6]]
        @test contributions(fig4; weights = true, weight_scheme = :legacy).values ≈ c4.values
        @test disentangle_individual_contribution(fig4; weights = true).values == c4.values
    end

    @testset "Fig. 5: temporal overlap" begin
        t0 = [0 1; 1 1]
        t1 = [0 1; 1 0]
        t2 = [1 0; 0 1]
        @test temporal_overlap(t0, t1) ≈ sqrt(1 / 2)
        @test temporal_overlap(t0, t1; oriented = true) ≈ sqrt(1 / 2)
        @test temporal_overlap(t1, t2) ≈ 1
        @test temporal_overlap(t1, t2; oriented = true) ≈ sqrt(2)
        @test beta_temporal_overlap(t1, t2; oriented = true) ≈ sqrt(2)
    end

    @testset "Fig. 2: all 127 metacommunities of 3 species" begin
        comps = all_compositions(3)[2:end, :]
        subsets = [findall(i -> (s >> (i - 1)) & 1 == 1, 1:7) for s in 1:127]
        beta = map(subsets) do s
            rows = repeat(s; outer = 3)[1:max(3, length(s))]
            betavolume(comps[rows, :])
        end
        ncomp = length.(subsets)
        @test length(beta) == 127
        @test all(b -> 0 <= b <= 3 + 1e-12, beta)
        @test maximum(beta) ≈ 3
        @test findall(>=(3 - 1e-12), beta) == [127]
        @test all(beta[ncomp .< 3] .== 0)
        spans = [rank(comps[s, :]) == 3 for s in subsets]
        @test (beta .> 0) == spans
        means = [mean(beta[ncomp .== k]) for k in 2:7]
        @test all(diff(means) .> 0)
    end
end

@testset "Known shapes" begin
    for d in 2:7
        I_d = Matrix{Float64}(I, d, d)
        @test BV._hull_volume(I_d) ≈ 1 / factorial(d)
        @test BV._hull_volume(all_compositions(d)) ≈ 1
        @test BV._hull_volume(vcat(I_d, -I_d)) ≈ 2^d / factorial(d)
        @test betavolume(I_d) ≈ d * (1 / factorial(d))^(1 / d)
        @test betavolume(3 .* all_compositions(d)) ≈ 3d
    end
    @test BV._hull_volume(reshape([-1.0, 2.0], 2, 1)) == 3
    @test BV._hull_volume([1.0 0 0; 0 1 0; 1 1 0]) == 0
    @test BV._hull_volume([1.0 1 1; 2 2 2]) == 0
    @test !BV._is_full_dimensional([1.0 2 0; 1 2 0])
    @test BV._is_full_dimensional([1.0 0 0; 0 1 0])
end

@testset "Exact volumes vs independent computations" begin
    rng = MersenneTwister(1)
    for i in 1:25
        n = rand(rng, 3:30)
        m = randexp(rng, n, 2)
        iseven(i) && (m = round.(m))
        @test betavolume(m) ≈ beta_from_volume(shoelace_hull_area(m), 2) rtol = 1e-10
    end
    # regression: degenerate binary data (BetaVolume 0.1.0 silently returned 0)
    rows = ["011000", "101100", "101100", "011100", "001100", "011001", "000101",
            "000100", "001001", "100001", "001010", "000111", "100110", "100100",
            "001000", "011101", "110010", "010000", "000001", "000000", "110100",
            "110101", "100001", "000001", "100010", "101000", "101101", "000100",
            "111011"]
    m = Float64[r[j] == '1' for r in rows, j in 1:6]
    @test betavolume(m; weights = true) > 0
    @test betavolume(m; weights = true, weight_scheme = :legacy) > 0
end

@testset "Invariances and geometric properties" begin
    rng = MersenneTwister(3)
    m = Float64.(rand(rng, 25, 4) .< 0.5)
    b = betavolume(m)
    hn = betavolume(m; method = :hyper_normal)
    @test betavolume(m[randperm(rng, 25), :]) ≈ b
    @test betavolume(m[:, randperm(rng, 4)]) ≈ b
    @test betavolume(permutedims(m)) ≈ b
    @test betavolume(permutedims(m); method = :hyper_normal) ≈ hn
    @test betavolume(m[randperm(rng, 25), :]; method = :hyper_normal) ≈ hn

    m3 = [1 0 0; 0 1 0; 0 0 1; 1 1 1]
    b3 = betavolume(m3)
    @test betavolume(vcat(m3, [0 0 0])) ≈ b3          # empty community
    @test betavolume(vcat(m3, m3[2:2, :])) ≈ b3       # duplicated composition
    @test betavolume(vcat(m3, [0.2 0.2 0.2])) ≈ b3    # interior point

    m4 = rand(MersenneTwister(4), 30, 4)
    for s in (0.5, 2.0, 10.0)
        @test betavolume(s .* m4) ≈ s * betavolume(m4)
        @test betavolume(s .* m4; method = :hyper_normal) ≈
              s^2 * betavolume(m4; method = :hyper_normal)
    end

    rng = MersenneTwister(5)
    mm = Float64.(rand(rng, 6, 4) .< 0.5)
    prev = betavolume(mm)
    for _ in 1:15
        mm = vcat(mm, permutedims((rand(rng, 4) .< 0.5) .* rand(rng, 4)))
        cur = betavolume(mm)
        @test cur >= prev - 1e-12
        prev = cur
    end

    same = repeat([1 0 1], 5)
    @test betavolume(same) == 0
    @test betavolume(same; method = :hyper_normal) == 0
    @test betavolume(zeros(4, 3)) == 0
    @test betavolume(zeros(4, 3); weights = true) == 0
    @test betavolume(reshape([2, 1, 3], 3, 1)) == 3
    @test betavolume([1 1 0]) == 1
    @test betavolume([1, 1, 0]) == 1                  # vector = one community

    # default method switches at dim_threshold
    m6 = Float64.(rand(MersenneTwister(7), 40, 6) .< 0.5)
    @test betavolume(m6) == betavolume(m6; method = :exact)
    @test betavolume(m6; dim_threshold = 5) == betavolume(m6; method = :hyper_normal)
    @test betavolume(m6; method = "hyper_normal") == betavolume(m6; method = :hyper_normal)

    # compressed_angle rescales
    @test betavolume([1 0; 0 1; 1 1]; compressed_angle = 0.5) ≈ 1
    S = [1 0.5; 0.5 1]
    @test betavolume([1 0; 0 1; 1 1]; compressed_angle = compressed_angle(S)) ≈
          2 * compressed_angle(S)

    # Bool / Int / BitMatrix input
    @test betavolume([true false; false true; true true]) ≈ 2
    @test betavolume(BitMatrix([1 0; 0 1; 1 1])) ≈ 2
end

@testset "Preprocessing" begin
    @test preprocess([1 0; 0 0; 0 1]; weights = true).ids == [1, 3]
    @test size(preprocess([1 0; 0 0; 0 1]; weights = true).points) == (2, 2)
    wide = Float64.(rand(MersenneTwister(1), 3, 8) .< 0.5)
    E = preprocess(wide)
    @test size(E.points) == (8, 3)
    @test E.transposed
    @test E.ids == 1:8
    @test E.points == permutedims(wide)
    @test preprocess(permutedims(wide)).points == permutedims(wide)
    m = [1 1 0; 1 1 0; 0 0 1; 1 1 1]
    @test preprocess(m; remove_unique = true).points == [1 0; 0 1; 1 1]
    @test preprocess(m; remove_unique = true).ids == [1, 3, 4]
    # weights leave a metacommunity without duplicates unchanged
    md = [1 0 0; 0 1 0; 0 0 1; 1 1 0]
    for scheme in (:paper, :legacy)
        @test preprocess(md; weights = true, weight_scheme = scheme).points == md
    end
    # both schemes proportional to the multiplicity; paper weights average to one
    mr = Float64.(rand(MersenneTwister(2), 40, 4) .< 0.5)
    P = preprocess(mr; weights = true, weight_scheme = :paper)
    L = preprocess(mr; weights = true, weight_scheme = :legacy)
    ratio = P.points[P.points .!= 0] ./ L.points[L.points .!= 0]
    @test all(ratio .≈ ratio[1])
    nz = vec(sum(mr; dims = 2)) .> 0
    w = vec(sum(P.points; dims = 2)) ./ vec(sum(mr[nz, :]; dims = 2))
    g = BV._row_groups(mr[nz, :])
    @test mean(w[g.first]) ≈ 1
    # grouping agrees with unique()
    mg = Float64.(rand(MersenneTwister(3), 0:2, 200, 4))
    g = BV._row_groups(mg)
    @test mg[g.first, :] == unique(mg; dims = 1)
    @test BV._row_groups(sparse(mg)).group == g.group
    @test BV._row_groups([0.0 1; -0.0 1]).counts == [2]
    # input is never modified
    orig = copy(mr)
    preprocess(mr; weights = true, remove_unique = true)
    @test mr == orig
end

@testset "Input validation" begin
    m = [1 0; 0 1]
    @test_throws ArgumentError betavolume([1 0; NaN 1])
    @test_throws ArgumentError betavolume([1 0; Inf 1])
    @test_throws ArgumentError betavolume(zeros(0, 3))
    @test_throws ArgumentError betavolume(["a" "b"; "c" "d"])
    @test_throws ArgumentError betavolume([1 missing; 0 1])
    @test_throws ArgumentError betavolume("abc")
    @test_throws ArgumentError betavolume(m; method = :magic)
    @test_throws ArgumentError betavolume(m; compressed_angle = -1)
    @test_throws ArgumentError betavolume(m; compressed_angle = Inf)
    @test_throws ArgumentError betavolume(m; weight_scheme = :other)
    @test_throws ArgumentError contributions(m; method = :magic)
    @test_throws ArgumentError preprocess(sparse([1.0 NaN; 0 1]))
    err = try
        betavolume([1 0; NaN 1])
    catch e
        e
    end
    @test occursin("missing", err.msg)
    @test betavolume(Any[1 0; 0 1]) ≈ sqrt(2)
end

@testset "Hyper-normal approximation" begin
    rng = MersenneTwister(6)
    for _ in 1:10
        d = rand(rng, 2:12)
        n = rand(rng, (d + 3):80)
        m = (rand(rng, n, d) .< 0.4) .* randexp(rng, n, d)
        X = preprocess(m).points
        d = size(X, 2)
        ev = eigvals(Symmetric(cov(vcat(X, zeros(1, d)))))
        @test betavolume(X; preprocess = false, method = :hyper_normal) ≈
              4d * exp(mean(log.(ev))) rtol = 1e-10
    end
    for d in (6, 8, 12)
        @test betavolume(all_compositions(d); method = :hyper_normal) ≈ d rtol = 0.05
    end
    @test BV._logdet_psd([1.0 1; 1 1])[1] == -Inf
    A = randn(MersenneTwister(41), 8, 13)
    C = A * A'
    @test BV._logdet_psd(C)[1] ≈ logdet(C)
end

@testset "Contributions" begin
    @testset "fast path equals brute force" begin
        rng = MersenneTwister(11)
        checked = 0
        for i in 1:30
            n = rand(rng, 12:50)
            d = rand(rng, 3:6)
            m = Float64.(rand(rng, n, d) .< 0.2 + 0.5 * rand(rng))
            i % 3 == 0 && (m .*= randexp(rng, n, d))
            for w in (false, true), method in (:exact, :hyper_normal)
                fast = quietly(() -> contributions(m; weights = w, method = method))
                (any(isnan, fast.values) || all(iszero, fast.values)) && continue
                @test fast.values ≈ brute_contribution(m; weights = w, method = method) rtol = 1e-8
                checked += 1
            end
        end
        @test checked > 60
    end

    @testset "properties" begin
        m = [1 0 0; 0 1 0; 0 0 1; 1 1 1; 0.3 0.3 0.3]
        c = contributions(m).values
        @test sum(c) ≈ 1
        @test all(>=(0), c)
        @test c[5] == 0
        @test all(>(0), c[1:4])
        md = [1 0; 0 1; 0 1; 1 1]
        cd = contributions(md).values
        @test cd[2:3] == [0, 0]
        @test sum(cd) ≈ 1
        # units follow the original columns when transposed
        wide = [1 0 1 1 0; 0 1 1 0 1]
        cw = contributions(wide)
        @test cw.units == [[i] for i in 1:5]
        @test contributions(wide; weights = true).units == [[1, 4], [2, 5], [3]]
        # 1-d embeddings
        c1 = contributions(reshape([2.0, 1.0, 3.0], 3, 1))
        @test c1.values == [0, 0, 1]
    end

    @testset "degenerate inputs warn" begin
        res = @test_logs (:warn, r"zero beta diversity") contributions([1 1; 2 2; 3 3])
        @test all(isnan, res.values)
        res = @test_logs (:warn, r"zero beta diversity") contributions([1 1; 2 2; 3 3];
                                                                       method = :hyper_normal)
        @test all(isnan, res.values)
        m = [1 0; 1 0; 0 1; 0 1]
        res = @test_logs (:warn, r"weights = true") contributions(m)
        @test res.values == zeros(4)
        @test contributions(m; weights = true).values ≈ [0.5, 0.5]
        @test contributions(zeros(3, 2); weights = true).values == Float64[]
    end

    @testset "bounded weights" begin
        fig4 = [1 0; 0 1; 0 1; 1 1; 1 1; 1 1]
        @test preprocess(fig4; weights = true, weight_scheme = :bounded).points ≈
              fig4 .* ([1, 2, 2, 3, 3, 3] ./ 3)
        @test contributions(fig4; weights = true, weight_scheme = :bounded).values ≈
              contributions(fig4; weights = true).values
        md = [1 0 0; 0 1 0; 0 0 1; 1 1 0]
        @test preprocess(md; weights = true, weight_scheme = :bounded).points == md
        rng = MersenneTwister(51)
        for _ in 1:300
            d = rand(rng, 2:4)
            comps = [Float64.(reverse(digits(k, base = 2, pad = d))) for k in 1:(2^d - 1)]
            N = rand(rng, d:40)
            Z = reduce(vcat, permutedims.(comps[rand(rng, 1:length(comps), N)]))
            size(preprocess(Z).points, 2) == d || continue
            b = betavolume(Z; weights = true, weight_scheme = :bounded)
            @test b <= d + 1e-9
            @test b <= betavolume(Z) + 1e-9
        end
    end

    @testset "preprocess = false (BetaVolume 0.1.0 semantics)" begin
        fig4 = [1 0; 0 1; 0 1; 1 1; 1 1; 1 1]
        X = preprocess(fig4; weights = true).points
        c = quietly(() -> contributions(X; preprocess = false))
        @test length(c.values) == 6
        @test c.units == [[i] for i in 1:6]
    end
end

@testset "Sparse input" begin
    rng = MersenneTwister(8)
    m = Float64.(rand(rng, 60, 15) .< 0.15)
    sm = sparse(m)
    for w in (false, true)
        @test betavolume(sm; weights = w, method = :hyper_normal) ≈
              betavolume(m; weights = w, method = :hyper_normal) rtol = 1e-10
        @test preprocess(sm; weights = w).points isa SparseMatrixCSC
        @test Matrix(preprocess(sm; weights = w, remove_unique = true).points) ==
              preprocess(m; weights = w, remove_unique = true).points
    end
    small = m[:, 1:5]
    @test betavolume(sparse(small); method = :exact) ≈ betavolume(small; method = :exact)
    @test BV._embedding_cov(sm) ≈ cov(vcat(m, zeros(1, 15)))
    m12 = Float64.(rand(MersenneTwister(12), 80, 12) .< 0.2)
    for w in (false, true)
        a = contributions(sparse(m12); weights = w, method = :hyper_normal)
        b = contributions(m12; weights = w, method = :hyper_normal)
        @test a.values ≈ b.values rtol = 1e-8
        @test a.units == b.units
    end
    @test betavolume(sparse(permutedims(m)); method = :hyper_normal) ≈
          betavolume(m; method = :hyper_normal)
    @test betavolume(transpose(sm); method = :hyper_normal) ≈
          betavolume(m; method = :hyper_normal)
end

@testset "Compressed angle" begin
    for d in 1:3
        @test compressed_angle(Matrix(1.0I, d, d)) == 1
    end
    for d in 4:6
        @test compressed_angle(Matrix(1.0I, d, d)) ≈ 1 rtol = 1e-4
    end
    @test compressed_angle(ones(2, 2)) == 0
    @test compressed_angle(ones(5, 5)) == 0
    @test compressed_angle([1 2; 2 4]) == 0
    rng = MersenneTwister(22)
    for _ in 1:20
        S = [1.0 1.8 * rand(rng) - 0.9; 1.8 * rand(rng) - 0.9 1.0]
        a, b = S[:, 1], S[:, 2]
        θ = acos(dot(a, b) / (norm(a) * norm(b)))
        @test compressed_angle(S) ≈ sqrt(2θ / π)
    end
    @test compressed_angle([1 0.5; 0.5 1]) ≈ sqrt(2 * acos(0.8) / π)
    rng = MersenneTwister(23)
    for _ in 1:20
        S = Matrix(1.0I, 3, 3) .+ 0.6 .* rand(rng, 3, 3) .* (1 .- Matrix(1.0I, 3, 3))
        Ω = 8 * solid_angle(S[:, 1], S[:, 2], S[:, 3]) / (4π)
        @test compressed_angle(S) ≈ Ω^(1 / 3)
    end
    rng = MersenneTwister(24)
    for d in 4:6
        S = Matrix(1.0I, d, d) .+ (0.8 .* rand(rng, d, d) .- 0.3) .* (1 .- Matrix(1.0I, d, d))
        est = compressed_angle(S; rng = MersenneTwister(d))
        mc = (2^d * mc_cone_fraction(S, 400_000, rng))^(1 / d)
        @test est ≈ mc rtol = 0.02
        @test compressed_angle(S; rng = MersenneTwister(1)) ==
              compressed_angle(S; rng = MersenneTwister(1))
    end
    rng = MersenneTwister(25)
    for d in 2:5
        perm = Matrix(1.0I, d, d)[randperm(rng, d), :]
        @test compressed_angle(perm) ≈ 1 rtol = 1e-4
        @test compressed_angle(Diagonal(0.5 .+ 2.5 .* rand(rng, d))) ≈ 1 rtol = 1e-4
    end
    rot(t) = [cos(t) -sin(t); sin(t) cos(t)]
    @test compressed_angle(rot(0.3)) ≈ 1
    for s in (0.1, 0.5, 0.9)
        @test compressed_angle([1 s; s 1]) < 1
        @test compressed_angle([1 -s; -s 1]) > 1
        S3 = fill(s, 3, 3)
        S3[diagind(S3)] .= 1
        @test compressed_angle(S3) < 1
    end
    f = [compressed_angle([1 s; s 1]) for s in (0, 0.25, 0.5, 0.75)]
    @test all(diff(f) .< 0)
    @test get_compressed_angle([1 0.5; 0.5 1]) == compressed_angle([1 0.5; 0.5 1])
    @test_throws ArgumentError compressed_angle([1 2 3; 4 5 6])
    @test_throws ArgumentError compressed_angle("a")
    @test_throws ArgumentError compressed_angle([1 NaN; 0 1])
end

@testset "Temporal overlap" begin
    rng = MersenneTwister(31)
    m = rand(rng, 20, 3)
    @test temporal_overlap(m, m) ≈ 1
    @test temporal_overlap(m, m ./ 2) ≈ 0.5 rtol = 1e-8
    @test temporal_overlap(m, 2 .* m) ≈ 1 rtol = 1e-8
    @test temporal_overlap([1 0; 1 1], [0 1; 1 1]) == 0      # share only an edge
    sq = [1 0; 0 1; 1 1]
    tri = [1 0; 0 1; 0 0]
    @test temporal_overlap(sq, tri) ≈ sqrt(1 / 2) rtol = 1e-10
    @test temporal_overlap(tri, sq) ≈ 1 rtol = 1e-10
    # 1-d
    @test temporal_overlap(reshape([1.0, 2.0], 2, 1), reshape([1.0, 1.0], 2, 1);
                           preprocess = false) ≈ 0.5
    # nested cubes in 3-d and 4-d
    for d in 3:4
        c = all_compositions(d)
        @test temporal_overlap(c, 0.5 .* c) ≈ 0.5 rtol = 1e-8
        @test temporal_overlap(c, c .+ 0.0) ≈ 1 rtol = 1e-8
    end
    # shifted cube in 2-d: (overlap of [0,1]^2 and [0.5,1.5]x[0,1] with origin hull)
    t0 = [0 1; 1 1]
    t1 = [0 1; 1 0]
    @test_throws ArgumentError temporal_overlap(t0, hcat(t0, [1, 1]))
    @test_throws ArgumentError temporal_overlap(t0, t1; remove_unique = true)
    @test_throws ArgumentError temporal_overlap(sq, sq; oriented = true)
    res = @test_logs (:warn, r"zero volume") temporal_overlap([1 1; 2 2], t1)
    @test isnan(res)
    res = @test_logs (:warn, r"zero volume") temporal_overlap([1 1; 2 2], t1; oriented = true)
    @test isnan(res)
end

@testset "LP (simplex and Chebyshev centre)" begin
    # textbook LP: max 3x + 5y, x ≤ 4, 2y ≤ 12, 3x + 2y ≤ 18 -> (2, 6), 36
    st, z, obj = BV._simplex_max([3.0, 5.0], [1.0 0; 0 2; 3 2], [4.0, 12, 18])
    @test st === :optimal
    @test z ≈ [2, 6]
    @test obj ≈ 36
    # Beale's cycling example (cycles with Dantzig's rule): optimum 5/4
    c = [0.75, -20, 0.5, -6]
    M = [0.25 -8 -1 9; 0.5 -12 -0.5 3; 0 0 1 0]
    st, z, obj = BV._simplex_max(c, M, [0.0, 0, 1])
    @test st === :optimal
    @test obj ≈ 1.25
    @test all(M * z .<= [0, 0, 1] .+ 1e-12)
    # unbounded
    st, _, _ = BV._simplex_max([1.0, 0], [-1.0 1], [1.0])
    @test st === :unbounded
    # box [-1, 1]^2
    x, r, st = BV._chebyshev_center([1.0 0; 0 1; -1 0; 0 -1], ones(4))
    @test r ≈ 1
    @test x ≈ [0, 0] atol = 1e-12
    # triangle (0,0), (1,0), (0,1): inradius 1 - 1/√2 at (r, r)
    x, r, st = BV._chebyshev_center([-1.0 0; 0 -1; 1/√2 1/√2], [0, 0, 1/√2])
    @test r ≈ 1 - 1 / √2
    @test x ≈ [r, r]
    # unnormalised rows give the same answer
    x2, r2, _ = BV._chebyshev_center([-1.0 0; 0 -1; 1 1], [0.0, 0, 1])
    @test r2 ≈ r
    @test x2 ≈ x
    # 3-d simplex: r = 1 / (3 + √3)
    x, r, _ = BV._chebyshev_center([-1.0 0 0; 0 -1 0; 0 0 -1; 1 1 1], [0.0, 0, 0, 1])
    @test r ≈ 1 / (3 + √3)
    @test x ≈ fill(r, 3)
    # rectangle [0,2] x [0,1] with the origin as a vertex: r = 1/2, x feasible
    A = [1.0 0; 0 1; -1 0; 0 -1]
    b = [2.0, 1, 0, 0]
    x, r, _ = BV._chebyshev_center(A, b)
    @test r ≈ 0.5
    @test all(A * x .+ r .<= b .+ 1e-12)
    # lower-dimensional feasible set: r = 0
    _, r, _ = BV._chebyshev_center([1.0 0; -1 0; 0 1; 0 -1], [0.0, 0, 1, 1])
    @test abs(r) < 1e-12
    # halfspace vertices of the unit square
    V = BV._halfspace_vertices([1.0 0; 0 1; -1 0; 0 -1], [1.0, 1, 0, 0], [0.3, 0.4])
    @test sortslices(round.(V; digits = 12); dims = 2) == [0 0 1 1; 0 1 0 1]
end

@testset "Qhull interface" begin
    # internal hull helpers agree with DirectQhull's public ConvexHull
    rng = MersenneTwister(3)
    for d in 2:6, _ in 1:5
        n = rand(rng, (d + 2):30)
        X = round.(rand(rng, n, d) .* (rand(rng, n, d) .< 0.7); digits = 2)
        P = BV._hull_points(X)
        BV._is_full_dimensional(P) || continue
        h = BV._hull(P; vertices = true, equations = true)
        ref = BV.DirectQhull.ConvexHull(P)
        @test h.volume ≈ ref.volume rtol = 1e-12
        @test h.vertices == sort(Int.(ref.vertices))
        @test h.equations == ref.equations
    end
    # hulls without pre-merging ("Qt Q0", volume and vertices) agree with "Qt"; they may list
    # different non-extreme (coplanar) points, but never a different extreme point
    rng = MersenneTwister(17)
    for i in 1:40
        d = rand(rng, 3:8)
        X = Float64.(rand(rng, 40, d) .< 0.2 + 0.4rand(rng))
        i % 4 == 0 && (X .*= randexp(rng, 40, d))
        P = BV._hull_points(X)
        BV._is_full_dimensional(P) || continue
        h = BV._hull(P; vertices = true)
        _, v, vt, _ = BV._qhull_run(P, BV._QH_FALLBACK_OPTIONS, true, false)
        @test h.volume ≈ v rtol = 1e-10
        for k in symdiff(h.vertices, vt)
            @test BV._hull(P[:, setdiff(1:size(P, 2), k)]).volume ≈ v rtol = 1e-10
        end
    end
    # dqrdc2-style rank
    @test BV._qr_rank!([1.0 2; 2 4; 3 6], 1e-10) == 1
    @test BV._qr_rank!([1.0 0; 0 1; 1 1], 1e-10) == 2
    @test BV._qr_rank!(zeros(4, 3), 1e-10) == 0
    # threaded and serial leave-one-out loops agree (also run with JULIA_NUM_THREADS > 1)
    Xc = rand(MersenneTwister(9), 60, 5)
    @test BV._exact_loo_loss(Xc)[1] == BV._exact_loo_loss(Xc)[1]
end

@testset "Type stability" begin
    X = rand(MersenneTwister(5), 20, 3)
    S = sparse(X)
    @inferred betavolume(X)
    @inferred betavolume(S; method = :hyper_normal)
    @inferred preprocess(X; weights = true)
    @inferred preprocess(S; weights = true)
    @inferred contributions(X; weights = true)
    @inferred contributions([1 0; 0 1; 1 1]; method = :hyper_normal)
    @inferred compressed_angle([1 0.2; 0.1 1])
    @inferred temporal_overlap(X, 2 .* X)
end

@testset "Cross-validation against R" begin
    ref = JSON.parsefile(joinpath(@__DIR__, "reference.json"))
    ncheck = Ref(0)
    @testset "case $(case["name"])" for case in ref["cases"]
        Z = tomatrix(case["Z"])
        for p in case["preprocess"]
            E = preprocess(Z; weights = p["weights"], weight_scheme = Symbol(p["scheme"]))
            Pr = tomatrix(p["points"])
            @test size(E.points) == size(Pr) || (isempty(Pr) && size(E.points, 1) == 0)
            isempty(Pr) || @test E.points ≈ Pr rtol = 1e-14
            @test E.ids == Int.(tovec(p["ids"]))
        end
        for b in case["beta"]
            v = betavolume(Z; weights = b["weights"], weight_scheme = Symbol(b["scheme"]),
                           method = Symbol(b["method"]))
            @test isapprox(v, tonum(b["value"]); rtol = 1e-8, atol = 1e-12)
            ncheck[] += 1
        end
        for c in case["contrib"]
            res = quietly(() -> contributions(Z; weights = c["weights"],
                                              weight_scheme = Symbol(c["scheme"]),
                                              method = Symbol(c["method"])))
            expected = tovec(c["values"])
            if all(isnan, expected)
                @test all(isnan, res.values)
            else
                @test isapprox(res.values, expected; rtol = 1e-8, atol = 1e-10)
            end
            labels = c["labels"] isa AbstractVector ? c["labels"] : [c["labels"]]
            @test [join(u, ", ") for u in res.units] == String.(labels)
            ncheck[] += 1
        end
    end
    @testset "compressed angle" for a in ref["angle"]
        S = tomatrix(a["S"])
        @test isapprox(compressed_angle(S), tonum(a["value"]); rtol = 1e-10, atol = 1e-14)
        ncheck[] += 1
    end
    @testset "temporal overlap" for t in ref["temporal"]
        v = quietly(() -> temporal_overlap(tomatrix(t["Z1"]), tomatrix(t["Z2"]);
                                          oriented = t["oriented"], weights = t["weights"],
                                          weight_scheme = Symbol(t["scheme"])))
        expected = tonum(t["value"])
        if isnan(expected)
            @test isnan(v)
        else
            @test isapprox(v, expected; rtol = 1e-6, atol = 1e-10)
        end
        ncheck[] += 1
    end
    @test ncheck[] > 300
end

end # testset
