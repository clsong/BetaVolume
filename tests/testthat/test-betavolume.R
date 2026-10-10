# Worked examples of Song et al. (2025) ----------------------------------------

test_that("Fig. 1: metacommunity I and II", {
  meta_I <- rbind(c(1, 0), c(0, 1))
  meta_II <- rbind(c(1, 0), c(0, 1), c(1, 1))
  # vol = 1/2 and 1  =>  beta = 2 * sqrt(vol)
  expect_equal(betavolume(meta_I), 2 * sqrt(1 / 2))
  expect_equal(betavolume(meta_II), 2)
  # Adding a distinct, species-rich community increases beta diversity
  expect_gt(betavolume(meta_II), betavolume(meta_I))
})

test_that("maximum beta diversity d is reached with all 2^d - 1 compositions", {
  for (d in 2:5) {
    full <- all_compositions(d)
    expect_equal(betavolume(full[-1, , drop = FALSE], hypervolume_method = "exact"), d)
    # the empty community does not matter
    expect_equal(betavolume(full, hypervolume_method = "exact"), d)
  }
})

test_that("Fig. 2: all metacommunities with three species", {
  comps <- all_compositions(3)[-1, ] # 7 non-empty compositions
  subsets <- unlist(lapply(1:7, function(k) combn(7, k, simplify = FALSE)),
    recursive = FALSE
  )
  # Communities are repeated so that N >= gamma and the embedding dimension is
  # gamma = 3 (duplicated compositions do not change the hull).
  beta <- vapply(subsets, function(s) {
    betavolume(comps[rep(s, length.out = max(3, length(s))), , drop = FALSE])
  }, numeric(1))
  n_comp <- lengths(subsets)

  expect_length(beta, 127)
  expect_true(all(beta >= 0 & beta <= 3 + 1e-12))
  # The maximum (3) is reached only with all seven compositions
  expect_equal(max(beta), 3)
  expect_identical(which(beta >= 3 - 1e-12), length(subsets))
  # Fewer than three compositions cannot span 3-d space
  expect_true(all(beta[n_comp < 3] == 0))
  # beta > 0 exactly when the compositions span R^3
  spans <- vapply(subsets, function(s) qr(comps[s, , drop = FALSE])$rank == 3, logical(1))
  expect_identical(beta > 0, spans)
  # Positive association with the number of unique compositions (Fig. 2)
  expect_true(all(diff(tapply(beta[n_comp >= 2], n_comp[n_comp >= 2], mean)) > 0))
})

# Exactness against independent computations ----------------------------------

test_that("2-d values match chull() + shoelace formula", {
  set.seed(1)
  for (i in 1:25) {
    n <- sample(3:30, 1)
    m <- matrix(stats::rexp(2 * n), n, 2)
    if (i %% 2 == 0) m <- round(m) # many duplicates / collinear points
    expected <- beta_from_volume(shoelace_hull_area(m), 2)
    expect_equal(betavolume(m), expected, tolerance = 1e-10)
  }
})

test_that("exact volumes match a Delaunay triangulation in 3-6 dimensions", {
  set.seed(2)
  for (i in 1:20) {
    d <- sample(3:6, 1)
    n <- sample((d + 2):40, 1)
    m <- matrix(stats::rbinom(n * d, 1, 0.5), n, d)
    if (i %% 2 == 0) m <- m * matrix(stats::runif(n * d), n, d)
    X <- preprocess_meta_composition(m)
    if (ncol(X) != d || qr(X)$rank < d) next
    expected <- beta_from_volume(delaunay_volume(X), d)
    expect_equal(betavolume(X, preprocess = FALSE, hypervolume_method = "exact"),
      expected,
      tolerance = 1e-9
    )
  }
})

test_that("degenerate binary data no longer silently returns 0 (regression)", {
  # BetaVolume 0.1.0 used Qhull's 'Tv' option, which aborts with an internal
  # error on this presence/absence matrix; the error was caught and 0 returned.
  rows <- c(
    "011000", "101100", "101100", "011100", "001100", "011001", "000101",
    "000100", "001001", "100001", "001010", "000111", "100110", "100100",
    "001000", "011101", "110010", "010000", "000001", "000000", "110100",
    "110101", "100001", "000001", "100010", "101000", "101101", "000100",
    "111011"
  )
  m <- t(vapply(strsplit(rows, ""), as.numeric, numeric(6)))
  for (scheme in c("paper", "legacy")) {
    X <- preprocess_meta_composition(m, weights = TRUE, weight_scheme = scheme)
    b <- betavolume(m, weights = TRUE, weight_scheme = scheme)
    expect_gt(b, 0)
    expect_equal(b, beta_from_volume(delaunay_volume(X), 6), tolerance = 1e-9)
  }
})

test_that("known shapes: simplex and scaled cube", {
  for (d in 2:6) {
    simplex <- diag(d) # vol = 1 / d!
    expect_equal(betavolume(simplex), d * (1 / factorial(d))^(1 / d))
    cube <- all_compositions(d) * 3 # vol = 3^d
    expect_equal(betavolume(cube), 3 * d)
  }
})

# Geometric properties ---------------------------------------------------------

test_that("beta diversity is invariant to permutations and orientation", {
  set.seed(3)
  m <- matrix(stats::rbinom(25 * 4, 1, 0.5), 25, 4)
  b <- betavolume(m)
  expect_equal(betavolume(m[sample(nrow(m)), ]), b)
  expect_equal(betavolume(m[, sample(ncol(m))]), b)
  expect_equal(betavolume(t(m)), b)
  hn <- betavolume(m, hypervolume_method = "hyper_normal")
  expect_equal(betavolume(t(m), hypervolume_method = "hyper_normal"), hn)
  expect_equal(betavolume(m[sample(nrow(m)), ], hypervolume_method = "hyper_normal"), hn)
})

test_that("empty communities, duplicates and interior points do not change beta_vol", {
  m <- rbind(c(1, 0, 0), c(0, 1, 0), c(0, 0, 1), c(1, 1, 1))
  b <- betavolume(m)
  expect_equal(betavolume(rbind(m, 0)), b) # empty community (double-zero asymmetry)
  expect_equal(betavolume(rbind(m, m[2, ])), b) # duplicated composition
  expect_equal(betavolume(rbind(m, c(0.2, 0.2, 0.2))), b) # point inside the hull
})

test_that("scaling: beta_vol is linear, the hyper-normal estimator quadratic", {
  set.seed(4)
  m <- matrix(stats::runif(30 * 4), 30, 4)
  for (s in c(0.5, 2, 10)) {
    expect_equal(betavolume(s * m), s * betavolume(m))
    # det(Cov)^(1/d) scales like a variance
    expect_equal(
      betavolume(s * m, hypervolume_method = "hyper_normal"),
      s^2 * betavolume(m, hypervolume_method = "hyper_normal")
    )
  }
})

test_that("adding a community never decreases exact beta diversity", {
  set.seed(5)
  m <- matrix(stats::rbinom(6 * 4, 1, 0.5), 6, 4)
  previous <- betavolume(m)
  for (i in 1:15) {
    m <- rbind(m, stats::rbinom(4, 1, 0.5) * stats::runif(4))
    current <- betavolume(m)
    expect_gte(current, previous - 1e-12)
    previous <- current
  }
})

test_that("identical communities have zero beta diversity", {
  m <- matrix(c(1, 0, 1), nrow = 5, ncol = 3, byrow = TRUE)
  expect_equal(betavolume(m), 0)
  expect_equal(betavolume(m, hypervolume_method = "hyper_normal"), 0)
  expect_equal(betavolume(matrix(0, 4, 3)), 0)
})

test_that("one-dimensional embeddings are handled", {
  expect_equal(betavolume(matrix(c(2, 1, 3), ncol = 1)), 3)
  expect_equal(betavolume(matrix(c(1, 1, 0), nrow = 1)), 1) # one community
})

# Hyper-normal approximation -----------------------------------------------------

test_that("hyper-normal equals 4 d det(cov)^(1/d) (formula of BetaVolume 0.1.0)", {
  set.seed(6)
  for (i in 1:10) {
    d <- sample(2:12, 1)
    n <- sample((d + 3):80, 1)
    m <- matrix(stats::rbinom(n * d, 1, 0.4) * stats::rexp(n * d), n, d)
    X <- preprocess_meta_composition(m)
    d <- ncol(X)
    ev <- eigen(stats::cov(rbind(X, 0)), symmetric = TRUE, only.values = TRUE)$values
    expected <- 4 * d * exp(mean(log(ev)))
    expect_equal(betavolume(X, preprocess = FALSE, hypervolume_method = "hyper_normal"),
      expected,
      tolerance = 1e-10
    )
  }
})

test_that("hyper-normal is calibrated: ~d when all compositions are present", {
  for (d in c(6, 8, 12)) {
    hn <- betavolume(all_compositions(d), hypervolume_method = "hyper_normal")
    expect_equal(hn, d, tolerance = 0.05)
  }
})

test_that("default method switches at dim_threshold", {
  set.seed(7)
  m <- matrix(stats::rbinom(40 * 6, 1, .5), 40, 6)
  expect_equal(betavolume(m), betavolume(m, hypervolume_method = "deterministic"))
  expect_equal(betavolume(m, dim_threshold = 5), betavolume(m, hypervolume_method = "hyper_normal"))
  expect_equal(betavolume(m, hypervolume_method = "exact"), betavolume(m, hypervolume_method = "deterministic"))
})

test_that("sparse and dense matrices give identical results", {
  set.seed(8)
  m <- matrix(stats::rbinom(60 * 15, 1, 0.15), 60, 15)
  sm <- Matrix::Matrix(m, sparse = TRUE)
  expect_s4_class(sm, "sparseMatrix")
  for (w in c(FALSE, TRUE)) {
    expect_equal(
      betavolume(sm, weights = w, hypervolume_method = "hyper_normal"),
      betavolume(m, weights = w, hypervolume_method = "hyper_normal")
    )
  }
  small <- m[, 1:5]
  expect_equal(
    betavolume(Matrix::Matrix(small, sparse = TRUE), hypervolume_method = "exact"),
    betavolume(small, hypervolume_method = "exact")
  )
  # backward-compatible `sparsity` argument
  expect_equal(
    betavolume(m, sparsity = TRUE, hypervolume_method = "hyper_normal"),
    betavolume(m, hypervolume_method = "hyper_normal")
  )
})

test_that("compressed_angle rescales the result", {
  m <- rbind(c(1, 0), c(0, 1), c(1, 1))
  expect_equal(betavolume(m, compressed_angle = 0.5), 1)
  S <- matrix(c(1, 0.5, 0.5, 1), 2)
  expect_equal(
    betavolume(m, compressed_angle = get_compressed_angle(S)),
    2 * get_compressed_angle(S)
  )
})

# Inputs -------------------------------------------------------------------------

test_that("data frames and logical matrices are accepted", {
  m <- rbind(c(1, 0), c(0, 1), c(1, 1))
  expect_equal(betavolume(as.data.frame(m)), 2)
  expect_equal(betavolume(m == 1), 2)
})

test_that("invalid inputs give informative errors", {
  m <- rbind(c(1, 0), c(0, 1))
  expect_error(betavolume(rbind(m, c(NA, 1))), "missing")
  expect_error(betavolume(rbind(m, c(Inf, 1))), "infinite")
  expect_error(betavolume(matrix(letters[1:4], 2)), "numeric")
  expect_error(betavolume(data.frame(a = 1:2, b = c("x", "y"))), "numeric")
  expect_error(betavolume(list(1, 2)), "matrix")
  expect_error(betavolume(m, hypervolume_method = "magic"))
  expect_error(betavolume(m, compressed_angle = -1), "compressed_angle")
  expect_error(betavolume(m, weights = NA), "weights")
  expect_error(betavolume(m, weight_scheme = "other"))
})

# Regression against BetaVolume 0.1.0 --------------------------------------------

test_that("results of BetaVolume 0.1.0 are reproduced", {
  # Values computed with the original implementation (commit 6b2ebce).
  ex <- example_meta_composition
  expect_equal(betavolume(ex), 4.47817790657455, tolerance = 1e-10)
  expect_equal(betavolume(ex, remove_unique = TRUE), 4.47817790657455, tolerance = 1e-10)
  expect_equal(betavolume(ex, hypervolume_method = "hyper_normal"), 4.51708592656972,
    tolerance = 1e-10
  )
  # The old weights are available as weight_scheme = "legacy"
  expect_equal(betavolume(ex, weights = TRUE, weight_scheme = "legacy"), 2.46193067726437,
    tolerance = 1e-10
  )
  expect_equal(
    betavolume(ex, weights = TRUE, weight_scheme = "legacy", hypervolume_method = "hyper_normal"),
    1.84983850777852,
    tolerance = 1e-10
  )
})
