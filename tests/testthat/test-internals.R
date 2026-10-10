test_that("hull volumes of standard polytopes", {
  hull_volume <- BetaVolume:::hull_volume
  for (d in 2:7) {
    expect_equal(hull_volume(diag(d)), 1 / factorial(d)) # standard simplex
    expect_equal(hull_volume(all_compositions(d)), 1) # unit cube
    expect_equal(hull_volume(rbind(diag(d), -diag(d))), 2^d / factorial(d)) # cross-polytope
  }
  expect_equal(hull_volume(matrix(c(-1, 2), ncol = 1)), 3) # 1-d interval
})

test_that("flat point clouds have zero volume", {
  hull_volume <- BetaVolume:::hull_volume
  expect_equal(hull_volume(rbind(c(1, 0, 0), c(0, 1, 0), c(1, 1, 0))), 0)
  expect_equal(hull_volume(rbind(c(1, 1, 1), c(2, 2, 2))), 0)
  expect_false(BetaVolume:::is_full_dimensional(rbind(c(1, 1), c(2, 2), c(0, 0))))
})

test_that("log-determinant agrees with base::determinant", {
  set.seed(41)
  for (i in 1:10) {
    d <- sample(2:15, 1)
    a <- matrix(stats::rnorm(d * (d + 5)), d)
    m <- tcrossprod(a)
    expect_equal(BetaVolume:::logdet_psd(m), as.numeric(determinant(m)$modulus))
  }
  expect_identical(BetaVolume:::logdet_psd(matrix(c(1, 1, 1, 1), 2)), -Inf)
})

test_that("sparse covariance equals stats::cov", {
  set.seed(42)
  m <- matrix(stats::rbinom(50 * 8, 1, 0.3), 50, 8)
  expected <- stats::cov(rbind(m, 0))
  dimnames(expected) <- NULL
  expect_equal(BetaVolume:::embedding_covariance(Matrix::Matrix(m, sparse = TRUE)), expected)
  expect_equal(BetaVolume:::embedding_covariance(m), expected)
})

test_that("the magrittr pipe is re-exported", {
  expect_equal(example_meta_composition %>% betavolume(), betavolume(example_meta_composition))
})
