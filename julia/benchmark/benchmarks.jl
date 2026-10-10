# Simple timing script (not part of the test suite).
#
#   julia --project=julia julia/benchmark/benchmarks.jl            # time Julia
#   julia --project=julia julia/benchmark/benchmarks.jl <dir>      # also export CSVs
#   Rscript julia/benchmark/benchmarks.R <dir>                     # time R on the same data
#
# Use JULIA_NUM_THREADS=k to parallelise the exact leave-one-out loop of contributions().

using BetaVolume
using Random
using Logging

function timeit(f; reps = 5, warmup = true)
    warmup && f()             # compile / warm up
    best = Inf
    for _ in 1:reps
        best = min(best, @elapsed f())
    end
    return best
end

function writecsv(path, M)
    open(path, "w") do io
        for i in axes(M, 1)
            println(io, join(M[i, :], ","))
        end
    end
end

rng = Xoshiro(2025)
bin200x8 = Float64.(rand(rng, 200, 8) .< 0.5)
bin2000x200 = Float64.(rand(rng, 2000, 200) .< 0.2)
bin2000x100 = Float64.(rand(rng, 2000, 100) .< 0.2)

if length(ARGS) >= 1
    mkpath(ARGS[1])
    writecsv(joinpath(ARGS[1], "bin200x8.csv"), Int.(bin200x8))
    writecsv(joinpath(ARGS[1], "bin2000x200.csv"), Int.(bin2000x200))
    writecsv(joinpath(ARGS[1], "bin2000x100.csv"), Int.(bin2000x100))
end

results = [
    ("betavolume exact, 200x8 binary, weights = true",
     () -> betavolume(bin200x8; weights = true, method = :exact)),
    ("betavolume hyper_normal, 2000x200 binary",
     () -> betavolume(bin2000x200; method = :hyper_normal)),
    ("contributions exact, 200x8 binary, weights = true",
     () -> contributions(bin200x8; weights = true, method = :exact)),
    ("contributions hyper_normal, 2000x100 binary",
     () -> contributions(bin2000x100; method = :hyper_normal)),
]

println("Julia ", VERSION, ", threads = ", Threads.nthreads())
with_logger(NullLogger()) do
    contributions(bin200x8[1:30, 1:4]; weights = true, method = :exact)  # compile
    for (name, f) in results
        slow = name == results[3][1]   # ~100+ hull recomputations in 8-d: one run only
        t = timeit(f; reps = slow ? 1 : 5, warmup = !slow)
        println(rpad(name, 55), lpad(string(round(t * 1000; digits = 2)), 12), " ms")
    end
end
