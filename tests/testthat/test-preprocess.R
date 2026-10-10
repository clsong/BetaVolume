fig4 <- rbind(
  c(1, 0),
  c(0, 1), c(0, 1),
  c(1, 1), c(1, 1), c(1, 1)
)

test_that("Fig. 4: weighted embedding of duplicated compositions (Eq. 5)", {
  w <- preprocess_meta_composition(fig4, weights = TRUE, weight_scheme = "paper")
  # (no. unique compositions) x (no. communities with the composition) / (no. communities)
  expected <- rbind(
    c(3 * 1 / 6, 0),
    c(0, 3 * 2 / 6), c(0, 3 * 2 / 6),
    c(3 * 3 / 6, 3 * 3 / 6), c(3 * 3 / 6, 3 * 3 / 6), c(3 * 3 / 6, 3 * 3 / 6)
  )
  expect_equal(unname(w), expected)
  expect_equal(betavolume(fig4, weights = TRUE, weight_scheme = "paper"), 2 * sqrt(1.125))
})

test_that("legacy weights (the default) match BetaVolume 0.1.0", {
  expect_identical(
    preprocess_meta_composition(fig4, weights = TRUE),
    preprocess_meta_composition(fig4, weights = TRUE, weight_scheme = "legacy")
  )
  w <- preprocess_meta_composition(fig4, weights = TRUE, weight_scheme = "legacy")
  # N n_k / sum_l n_l^2 with N = 6 and sum n_l^2 = 1 + 4 + 9 = 14
  expect_equal(unname(w[, 1]), 6 * c(1, 0, 0, 3, 3, 3) / 14)
  expect_equal(betavolume(fig4, weights = TRUE, weight_scheme = "legacy"), 1.818275, tolerance = 1e-6)
})

test_that("both weighting schemes are proportional to the multiplicity", {
  set.seed(1)
  m <- matrix(stats::rbinom(40 * 4, 1, 0.5), 40, 4)
  p <- preprocess_meta_composition(m, weights = TRUE, weight_scheme = "paper")
  l <- preprocess_meta_composition(m, weights = TRUE, weight_scheme = "legacy")
  ratio <- p[p != 0] / l[l != 0]
  expect_equal(ratio, rep(ratio[1], length(ratio)))
  # Paper weights average to one over the unique compositions
  keys <- apply(m[rowSums(m) > 0, ], 1, paste, collapse = "")
  w <- rowSums(p) / rowSums(m[rowSums(m) > 0, ])
  expect_equal(mean(w[!duplicated(keys)]), 1)
})

test_that("weights leave a metacommunity without duplicates unchanged", {
  m <- rbind(c(1, 0, 0), c(0, 1, 0), c(0, 0, 1), c(1, 1, 0))
  for (scheme in c("paper", "legacy")) {
    expect_equal(unname(preprocess_meta_composition(m, weights = TRUE, weight_scheme = scheme)), m)
  }
})

test_that("weighting drops empty communities", {
  m <- rbind(c(1, 0), c(0, 0), c(0, 1))
  w <- preprocess_meta_composition(m, weights = TRUE)
  expect_equal(nrow(w), 2)
  expect_equal(rownames(w), c("1", "3"))
})

test_that("orientation: rows become the embedded points", {
  wide <- matrix(stats::rbinom(3 * 8, 1, 0.5), 3, 8, dimnames = list(paste0("site", 1:3), paste0("sp", 1:8)))
  p <- preprocess_meta_composition(wide)
  expect_equal(dim(p), c(8, 3))
  expect_equal(rownames(p), paste0("sp", 1:8))
  tall <- t(wide)
  expect_equal(preprocess_meta_composition(tall), tall)
  # unnamed rows receive their index as name
  expect_equal(rownames(preprocess_meta_composition(matrix(1:6, 3))), c("1", "2", "3"))
})

test_that("remove_unique keeps one copy of duplicated rows and columns", {
  m <- cbind(c(1, 1, 0, 1), c(1, 1, 0, 1), c(0, 0, 1, 1))
  p <- preprocess_meta_composition(m, remove_unique = TRUE)
  expect_equal(unname(p), rbind(c(1, 0), c(0, 1), c(1, 1)))
})

test_that("sparse input stays sparse and matches dense results", {
  set.seed(2)
  m <- matrix(stats::rbinom(30 * 5, 1, 0.3), 30, 5)
  sm <- Matrix::Matrix(m, sparse = TRUE)
  for (w in c(FALSE, TRUE)) for (u in c(FALSE, TRUE)) {
    ps <- preprocess_meta_composition(sm, weights = w, remove_unique = u)
    pd <- preprocess_meta_composition(m, weights = w, remove_unique = u)
    expect_s4_class(ps, "sparseMatrix")
    expect_equal(as.matrix(ps), pd)
  }
  expect_equal(dim(preprocess_meta_composition(Matrix::t(sm))), c(30, 5))
})

test_that("row keys agree with base::duplicated for dense and sparse matrices", {
  set.seed(3)
  m <- matrix(stats::rbinom(200 * 4, 2, 0.3), 200, 4)
  expect_identical(duplicated(BetaVolume:::row_keys(m)), as.vector(duplicated(m)))
  sm <- Matrix::Matrix(m, sparse = TRUE)
  expect_identical(duplicated(BetaVolume:::row_keys(sm)), as.vector(duplicated(m)))
})

test_that("legacy weights keep beta_vol within [0, d] more often than Eq. 5", {
  # Fig. 4: Eq. 5 weights exceed the maximum d = 2, legacy weights do not
  expect_gt(betavolume(fig4, weights = TRUE, weight_scheme = "paper"), 2)
  expect_lte(betavolume(fig4, weights = TRUE), 2)
})
