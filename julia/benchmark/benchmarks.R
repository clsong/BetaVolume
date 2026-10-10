# R timings on the matrices exported by benchmarks.jl (see that file).
#   Rscript julia/benchmark/benchmarks.R <dir>
args <- commandArgs(trailingOnly = TRUE)
dir <- if (length(args)) args[1] else "."
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
here <- dirname(normalizePath(sub("^--file=", "", file_arg[1])))
suppressMessages(devtools::load_all(normalizePath(file.path(here, "..", "..")), quiet = TRUE))

rd <- function(f) unname(as.matrix(read.csv(file.path(dir, f), header = FALSE)))
bin200x8 <- rd("bin200x8.csv")
bin2000x200 <- rd("bin2000x200.csv")
bin2000x100 <- rd("bin2000x100.csv")

timeit <- function(expr, reps = 5, warmup = TRUE) {
  f <- function() suppressWarnings(eval(expr, envir = globalenv()))
  if (warmup) f()
  min(vapply(seq_len(reps), function(i) system.time(f())[["elapsed"]], numeric(1)))
}
cat("R", R.version$major, ".", R.version$minor, "\n", sep = "")
res <- list(
  "betavolume exact, 200x8 binary, weights = true" =
    quote(betavolume(bin200x8, weights = TRUE, hypervolume_method = "exact")),
  "betavolume hyper_normal, 2000x200 binary" =
    quote(betavolume(bin2000x200, hypervolume_method = "hyper_normal")),
  "contributions exact, 200x8 binary, weights = true" =
    quote(disentangle_individual_contribution(bin200x8, weights = TRUE, hypervolume_method = "exact")),
  "contributions hyper_normal, 2000x100 binary" =
    quote(disentangle_individual_contribution(bin2000x100, hypervolume_method = "hyper_normal"))
)
for (nm in names(res)) {
  slow <- grepl("contributions exact", nm)
  t <- timeit(res[[nm]], reps = if (slow) 1 else 5, warmup = !slow)
  cat(formatC(nm, width = -55), formatC(t * 1000, format = "f", digits = 2, width = 12), " ms\n")
}
cat("values:", betavolume(bin200x8, weights = TRUE), betavolume(bin2000x200, hypervolume_method = "hyper_normal"), "\n")
