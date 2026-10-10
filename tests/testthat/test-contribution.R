fig4 <- rbind(c(1, 0), c(0, 1), c(0, 1), c(1, 1), c(1, 1), c(1, 1))
fig6 <- rbind(c(0.5, 0), c(0, 1), c(1.5, 1.5))

test_that("Fig. 6: contributions 0.16, 0.37 and 0.47", {
  # Hand calculation: vol(P0) = 1.125, vol(P1) = 0.75, vol(P2) = 0.375, vol(P3) = 0.25
  loss <- sqrt(1.125) - sqrt(c(0.75, 0.375, 0.25))
  expected <- loss / sum(loss)
  from_matrix <- disentangle_individual_contribution(fig6)
  expect_equal(unname(from_matrix), expected)
  expect_equal(round(unname(from_matrix), 2), c(0.16, 0.37, 0.47))

  # Same result from the raw presence/absence data with the weighted embedding
  from_raw <- disentangle_individual_contribution(fig4, weights = TRUE)
  expect_equal(unname(from_raw), expected)
  expect_named(from_raw, c("1", "2, 3", "4, 5, 6"))
  # Both weighting schemes give the same contributions (scale invariance)
  expect_equal(disentangle_individual_contribution(fig4, weights = TRUE, weight_scheme = "legacy"), from_raw)
})

test_that("fast leave-one-out equals brute-force recomputation", {
  set.seed(11)
  checked <- 0
  for (i in 1:30) {
    n <- sample(12:50, 1)
    d <- sample(3:6, 1)
    m <- matrix(stats::rbinom(n * d, 1, stats::runif(1, 0.2, 0.7)), n, d)
    if (i %% 3 == 0) m <- m * matrix(stats::rexp(n * d), n, d)
    for (w in c(FALSE, TRUE)) {
      for (method in c("deterministic", "hyper_normal")) {
        fast <- suppressWarnings(
          disentangle_individual_contribution(m, weights = w, hypervolume_method = method)
        )
        if (anyNA(fast) || all(fast == 0)) next # degenerate cases tested below
        expect_equal(unname(fast), brute_contribution(m, w, method), tolerance = 1e-8)
        checked <- checked + 1
      }
    }
  }
  expect_gt(checked, 60)
})

test_that("exact contributions are non-negative, sum to one and vanish inside the hull", {
  m <- rbind(c(1, 0, 0), c(0, 1, 0), c(0, 0, 1), c(1, 1, 1), c(0.3, 0.3, 0.3))
  ctb <- disentangle_individual_contribution(m)
  expect_equal(sum(ctb), 1)
  expect_true(all(ctb >= 0))
  expect_equal(unname(ctb[5]), 0) # interior point
  # unique compositions on the hull all contribute positively
  expect_true(all(ctb[1:4] > 0))
})

test_that("duplicated compositions contribute 0 without weights", {
  m <- rbind(c(1, 0), c(0, 1), c(0, 1), c(1, 1))
  ctb <- disentangle_individual_contribution(m)
  expect_equal(unname(ctb[2:3]), c(0, 0))
  expect_equal(sum(ctb), 1)
})

test_that("names follow row (or transposed column) names", {
  ctb <- disentangle_individual_contribution(example_meta_composition)
  expect_named(ctb, rownames(example_meta_composition))
  m <- rbind(a = c(1, 0), b = c(0, 1), c = c(1, 1))
  expect_named(disentangle_individual_contribution(m), c("a", "b", "c"))
})

test_that("degenerate inputs warn instead of returning garbage", {
  expect_warning(
    res <- disentangle_individual_contribution(rbind(c(1, 1), c(2, 2), c(3, 3))),
    "zero beta diversity"
  )
  expect_true(all(is.na(res)))
  # Every vertex occurs twice: no single row matters
  m <- rbind(c(1, 0), c(1, 0), c(0, 1), c(0, 1))
  expect_warning(res <- disentangle_individual_contribution(m), "weights = TRUE")
  expect_equal(unname(res), rep(0, 4))
  expect_equal(unname(disentangle_individual_contribution(m, weights = TRUE)), c(0.5, 0.5))
})

test_that("hyper-normal contributions work for sparse input", {
  set.seed(12)
  m <- matrix(stats::rbinom(80 * 12, 1, 0.2), 80, 12)
  sm <- Matrix::Matrix(m, sparse = TRUE)
  expect_equal(
    disentangle_individual_contribution(sm, hypervolume_method = "hyper_normal"),
    disentangle_individual_contribution(m, hypervolume_method = "hyper_normal")
  )
  expect_equal(
    disentangle_individual_contribution(sm, weights = TRUE, hypervolume_method = "hyper_normal"),
    disentangle_individual_contribution(m, weights = TRUE, hypervolume_method = "hyper_normal")
  )
})

test_that("results of BetaVolume 0.1.0 are reproduced", {
  # BetaVolume 0.1.0 removed the rows of the weighted embedding one at a time;
  # this corresponds to preprocess = FALSE on the weighted matrix.
  X <- preprocess_meta_composition(example_meta_composition, weights = TRUE)
  ctb <- disentangle_individual_contribution(X, preprocess = FALSE)
  top <- c(10, 21, 29, 22, 9, 12)
  expect_equal(unname(ctb[top]),
    c(0.211437843333, 0.139613964289, 0.122300871678, 0.111325198370, 0.093546806763, 0.090604770562),
    tolerance = 1e-9
  )
  expect_equal(sum(ctb == 0), 20)
})

test_that("parallel and serial exact contributions are identical", {
  skip_on_os("windows")
  skip_on_cran()
  set.seed(13)
  m <- matrix(stats::rbinom(40 * 5, 1, 0.5), 40, 5)
  expect_identical(
    disentangle_individual_contribution(m, weights = TRUE, cores = 2),
    disentangle_individual_contribution(m, weights = TRUE, cores = 1)
  )
  expect_error(disentangle_individual_contribution(m, cores = 0), "cores")
})

test_that("hyper-normal contributions do not overflow with many duplicates", {
  # n * n_k exceeded .Machine$integer.max in BetaVolume 0.2.0-dev
  m <- rbind(matrix(c(1, 0, 0), 46500, 3, byrow = TRUE), diag(3), c(1, 1, 0), c(0, 1, 1))
  ctb <- disentangle_individual_contribution(m,
    weights = TRUE, weight_scheme = "paper", hypervolume_method = "hyper_normal"
  )
  expect_false(anyNA(ctb))
  expect_equal(sum(ctb), 1)
})
