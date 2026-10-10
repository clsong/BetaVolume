# Benchmark BetaVolume 0.2.0 against BetaVolume 0.1.0 (commit 6b2ebce).
#
# Usage:
#   R CMD INSTALL -l <legacy_lib> <checkout of commit 6b2ebce>
#   Rscript benchmarks/benchmark_r.R <legacy_lib> [<repo>]
#
# run from the root of the repository. The two versions are loaded in
# separate Rscript processes and timed on identical inputs.

args <- commandArgs(trailingOnly = TRUE)
legacy_lib <- if (length(args) >= 1) args[[1]] else "~/R/legacy"
repo <- normalizePath(if (length(args) >= 2) args[[2]] else getwd())

cases <- quote({
  set.seed(1)
  bin_200x8 <- matrix(rbinom(200 * 8, 1, 0.4), 200, 8)
  bin_60x6 <- matrix(rbinom(60 * 6, 1, 0.4), 60, 6)
  abund_2000x200 <- matrix(rbinom(2000 * 200, 1, 0.1) * rexp(2000 * 200), 2000, 200)
  abund_1000x50 <- matrix(rbinom(1000 * 50, 1, 0.2) * rexp(1000 * 50), 1000, 50)
  ex <- example_meta_composition
  timeit <- function(expr, reps = 3) {
    expr <- substitute(expr)
    t <- vapply(seq_len(reps), function(i) {
      system.time(suppressWarnings(eval(expr, parent.frame(3))))[["elapsed"]]
    }, numeric(1))
    stats::median(t)
  }
  res <- c(
    "betavolume, example data, weights (exact)" = timeit(betavolume(ex, weights = TRUE), 20),
    "betavolume, 200x8 binary, weights (exact)" = timeit(betavolume(bin_200x8, weights = TRUE)),
    "betavolume, 2000x200, hyper_normal" = timeit(betavolume(abund_2000x200, hypervolume_method = "hyper_normal")),
    "contributions, example data, weights (exact)" = timeit(disentangle_individual_contribution(ex, weights = TRUE)),
    "contributions, 60x6 binary (exact)" = timeit(disentangle_individual_contribution(bin_60x6)),
    "contributions, 200x8 binary, weights (exact)" = timeit(disentangle_individual_contribution(bin_200x8, weights = TRUE), 1),
    "contributions, 1000x50 (hyper_normal)" = timeit(disentangle_individual_contribution(abund_1000x50), 1)
  )
  res
})

run_version <- function(lib_expr) {
  script <- tempfile(fileext = ".R")
  out <- tempfile(fileext = ".rds")
  writeLines(c(lib_expr, paste0("res <- ", paste(deparse(cases), collapse = "\n")), sprintf("saveRDS(res, '%s')", out)), script)
  log <- tempfile(fileext = ".log")
  system2("Rscript", script, stdout = log, stderr = log)
  if (!file.exists(out)) stop("benchmark run failed:\n", paste(readLines(log), collapse = "\n"))
  readRDS(out)
}

new <- run_version(sprintf("suppressMessages(devtools::load_all('%s', quiet = TRUE))", repo))
old <- run_version(sprintf("suppressWarnings(library(BetaVolume, lib.loc = '%s'))", legacy_lib))

tab <- data.frame(
  task = names(new),
  `v0.1.0 (s)` = signif(old, 3),
  `v0.2.0 (s)` = signif(new, 3),
  speedup = paste0(round(old / new, 1), "x"),
  check.names = FALSE, row.names = NULL
)
print(tab, right = FALSE)
if (requireNamespace("knitr", quietly = TRUE)) print(knitr::kable(tab))
