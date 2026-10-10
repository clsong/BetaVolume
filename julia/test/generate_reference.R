# Generate reference values from the R implementation for the Julia test suite.
#
# Usage (from anywhere):  Rscript julia/test/generate_reference.R
# Writes julia/test/reference.json next to this script. Requires devtools and jsonlite;
# the R package is loaded from the repository with devtools::load_all().

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
here <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
repo <- normalizePath(file.path(here, "..", ".."))
suppressMessages(devtools::load_all(repo, quiet = TRUE))

quiet <- function(expr) suppressWarnings(expr)

# Matrix generators (entries rounded so that they are exactly representable in JSON)
gen_binary <- function(n, p, prob) matrix(stats::rbinom(n * p, 1, prob), n, p)
gen_abund <- function(n, p, prob) {
  matrix(stats::rbinom(n * p, 1, prob) * round(stats::rexp(n * p, 0.5), 3), n, p)
}
gen_dups <- function(n, p, prob, pool) {
  base <- gen_binary(pool, p, prob)
  base[sample(pool, n, replace = TRUE), , drop = FALSE]
}
gen_counts <- function(n, p, lambda) matrix(stats::rpois(n * p, lambda), n, p)

set.seed(20250305)
specs <- list(
  list("binary", 12, 2, 0.5), list("abund", 15, 2, 0.7), list("binary", 20, 3, 0.4),
  list("abund", 25, 3, 0.6), list("dups", 30, 3, 0.5, 6), list("counts", 18, 3, 1.5),
  list("binary", 30, 4, 0.5), list("abund", 40, 4, 0.5), list("dups", 40, 4, 0.5, 9),
  list("binary", 35, 5, 0.4), list("abund", 30, 5, 0.6), list("dups", 50, 5, 0.5, 14),
  list("binary", 45, 6, 0.5), list("abund", 50, 6, 0.4), list("dups", 60, 6, 0.5, 20),
  list("counts", 40, 6, 1.0), list("binary", 50, 7, 0.5), list("abund", 45, 7, 0.5),
  list("binary", 60, 8, 0.4), list("abund", 40, 8, 0.6), list("binary", 70, 9, 0.5),
  list("abund", 60, 10, 0.4), list("binary", 80, 11, 0.3), list("abund", 80, 12, 0.5),
  list("dups", 80, 12, 0.3, 30),
  # wide matrices (transposed by preprocessing)
  list("binary", 5, 12, 0.5), list("abund", 4, 20, 0.6), list("binary", 8, 40, 0.3),
  list("abund", 3, 7, 0.7),
  # an empty row and an empty column
  list("binary_empty", 25, 4, 0.5), list("abund_empty", 20, 5, 0.5),
  # degenerate: rank deficient
  list("collinear", 10, 3, 0)
)

make_matrix <- function(s) {
  type <- s[[1]]; n <- s[[2]]; p <- s[[3]]; prob <- s[[4]]
  switch(type,
    binary = gen_binary(n, p, prob),
    abund = gen_abund(n, p, prob),
    dups = gen_dups(n, p, prob, s[[5]]),
    counts = gen_counts(n, p, prob),
    binary_empty = { m <- gen_binary(n, p, prob); m[3, ] <- 0; m[, 2] <- 0; m },
    abund_empty = { m <- gen_abund(n, p, prob); m[1, ] <- 0; m[, p] <- 0; m },
    collinear = outer(seq_len(n), c(1, 2, 3)) / n
  )
}

cases <- list()
for (i in seq_along(specs)) {
  m <- make_matrix(specs[[i]])
  d <- min(dim(m))
  betas <- list()
  contribs <- list()
  pre <- list()
  for (w in c(FALSE, TRUE)) {
    schemes <- if (w) c("paper", "legacy") else "paper"
    for (scheme in schemes) {
      X <- preprocess_meta_composition(m, weights = w, weight_scheme = scheme)
      pre[[length(pre) + 1]] <- list(
        weights = w, scheme = scheme, points = unname(as.matrix(X)),
        ids = as.integer(rownames(X))
      )
      methods <- if (d <= 8) c("exact", "hyper_normal") else "hyper_normal"
      for (method in methods) {
        b <- quiet(betavolume(m, weights = w, weight_scheme = scheme, hypervolume_method = method))
        betas[[length(betas) + 1]] <- list(method = method, weights = w, scheme = scheme, value = b)
      }
      cmethods <- if (d <= 6) c("exact", "hyper_normal") else "hyper_normal"
      for (method in cmethods) {
        ctb <- quiet(disentangle_individual_contribution(m,
          weights = w, weight_scheme = scheme, hypervolume_method = method
        ))
        contribs[[length(contribs) + 1]] <- list(
          method = method, weights = w, scheme = scheme,
          values = I(unname(as.numeric(ctb))), labels = I(names(ctb))
        )
      }
    }
  }
  cases[[i]] <- list(
    name = paste(c(specs[[i]][[1]], dim(m)), collapse = "_"),
    Z = m, d = d, beta = betas, contrib = contribs, preprocess = pre
  )
  message(sprintf("case %2d: %s (d = %d)", i, cases[[i]]$name, d))
}

# Compressed angle (closed forms, d <= 3)
angle_mats <- list(
  diag(2), matrix(c(1, 0.5, 0.5, 1), 2), matrix(c(1, -0.5, -0.5, 1), 2),
  matrix(c(1, 0.3, -0.2, 1), 2), matrix(c(1, 2, 2, 4), 2), matrix(2.5, 1, 1),
  matrix(c(1, 0.2, 0.4, 0.1, 1, 0.3, 0.5, 0.25, 1), 3),
  matrix(c(1, -0.3, 0.1, 0.2, 1, -0.4, 0.05, 0.3, 1), 3),
  matrix(c(1, 0.9, 0.9, 0.9, 1, 0.9, 0.9, 0.9, 1), 3), matrix(1, 3, 3)
)
angles <- lapply(angle_mats, function(S) {
  list(S = S, value = as.numeric(get_compressed_angle(S)))
})

# Temporal overlaps (2-4 dimensions)
temporal <- list()
add_temporal <- function(m1, m2, oriented = FALSE, weights = FALSE) {
  v <- quiet(beta_temporal_overlap(m1, m2, oriented = oriented, weights = weights))
  temporal[[length(temporal) + 1]] <<- list(
    Z1 = m1, Z2 = m2, oriented = oriented, weights = weights, value = v
  )
}
set.seed(77)
for (p in 2:4) {
  for (k in 1:3) {
    n <- sample((p + 2):12, 1)
    m1 <- gen_abund(n, p, 0.7)
    m2 <- if (k == 1) m1 * round(stats::runif(1, 0.5, 1.5), 2) else gen_abund(n, p, 0.7)
    add_temporal(m1, m2)
    b1 <- gen_binary(n, p, 0.5)
    b2 <- gen_binary(n, p, 0.5)
    add_temporal(b1, b2)
    add_temporal(b1, b2, weights = TRUE)
  }
  for (k in 1:3) {
    s1 <- matrix(round(stats::runif(p * p), 3), p, p)
    s2 <- if (k == 3) s1[, c(2, 1, seq_len(p)[-(1:2)])] else matrix(round(stats::runif(p * p), 3), p, p)
    add_temporal(s1, s2)
    add_temporal(s1, s2, oriented = TRUE)
  }
}

ref <- list(
  generated_with = paste("BetaVolume (R)", as.character(utils::packageVersion("BetaVolume")),
    "on", R.version.string),
  cases = cases, angle = angles, temporal = temporal
)
out <- file.path(here, "reference.json")
jsonlite::write_json(ref, out, digits = I(17), auto_unbox = TRUE, na = "null", pretty = FALSE)
message("wrote ", out)
