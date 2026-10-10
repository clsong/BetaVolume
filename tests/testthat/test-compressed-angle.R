test_that("identity similarity leaves the hypervolume unchanged", {
  for (d in 1:3) expect_identical(get_compressed_angle(diag(d)), 1)
  set.seed(21)
  for (d in 4:7) expect_equal(as.numeric(get_compressed_angle(diag(d))), 1, tolerance = 1e-4)
})

test_that("totally similar species or singular matrices collapse the hypervolume", {
  expect_identical(get_compressed_angle(matrix(1, 2, 2)), 0)
  expect_identical(get_compressed_angle(matrix(1, 5, 5)), 0)
  expect_identical(get_compressed_angle(rbind(c(1, 2), c(2, 4))), 0)
})

test_that("2-d values equal the planar angle between the transformed axes", {
  # The cone spanned by two columns covers theta / (2 pi) of the plane, so
  # Omega = 4 theta / (2 pi) and the factor is Omega^(1/2).
  set.seed(22)
  for (i in 1:20) {
    S <- diag(2)
    S[1, 2] <- stats::runif(1, -0.9, 0.9)
    S[2, 1] <- stats::runif(1, -0.9, 0.9)
    a <- S[, 1]
    b <- S[, 2]
    theta <- acos(sum(a * b) / sqrt(sum(a^2) * sum(b^2)))
    expect_equal(get_compressed_angle(S), sqrt(2 * theta / pi))
  }
  # Fig. 7C
  theta <- acos(0.8)
  expect_equal(get_compressed_angle(matrix(c(1, .5, .5, 1), 2)), sqrt(2 * theta / pi))
})

test_that("3-d values equal the Van Oosterom-Strackee solid angle", {
  solid_angle <- function(a, b, c) {
    num <- abs(sum(a * c(b[2] * c[3] - b[3] * c[2], b[3] * c[1] - b[1] * c[3], b[1] * c[2] - b[2] * c[1])))
    la <- sqrt(sum(a^2))
    lb <- sqrt(sum(b^2))
    lc <- sqrt(sum(c^2))
    den <- la * lb * lc + sum(a * b) * lc + sum(a * c) * lb + sum(b * c) * la
    om <- 2 * atan2(num, den)
    if (om < 0) om <- om + 4 * pi
    om
  }
  set.seed(23)
  for (i in 1:20) {
    S <- diag(3) + matrix(stats::runif(9, 0, 0.6), 3) * (1 - diag(3))
    omega <- 8 * solid_angle(S[, 1], S[, 2], S[, 3]) / (4 * pi)
    expect_equal(get_compressed_angle(S), omega^(1 / 3))
  }
})

test_that("higher dimensions agree with a Monte Carlo estimate", {
  set.seed(24)
  for (d in 4:6) {
    S <- diag(d) + matrix(stats::runif(d * d, -0.3, 0.5), d) * (1 - diag(d))
    est <- as.numeric(get_compressed_angle(S))
    mc <- (2^d * mc_cone_fraction(S, 4e5))^(1 / d)
    expect_equal(est, mc, tolerance = 0.02)
    expect_true(is.numeric(attr(get_compressed_angle(S), "error")))
  }
})

test_that("rotations, permutations and positive rescaling of axes give 1", {
  set.seed(25)
  for (d in 2:5) {
    perm <- diag(d)[sample(d), ]
    expect_equal(as.numeric(get_compressed_angle(perm)), 1, tolerance = 1e-4)
    scaled <- diag(stats::runif(d, 0.5, 3))
    expect_equal(as.numeric(get_compressed_angle(scaled)), 1, tolerance = 1e-4)
  }
  rot <- function(t) matrix(c(cos(t), sin(t), -sin(t), cos(t)), 2)
  expect_equal(get_compressed_angle(rot(0.3)), 1)
})

test_that("similarity shrinks and complementarity expands the hypervolume", {
  for (s in c(0.1, 0.5, 0.9)) {
    expect_lt(get_compressed_angle(matrix(c(1, s, s, 1), 2)), 1)
    expect_gt(get_compressed_angle(matrix(c(1, -s, -s, 1), 2)), 1)
    S3 <- matrix(s, 3, 3)
    diag(S3) <- 1
    expect_lt(get_compressed_angle(S3), 1)
  }
  # More similarity, more shrinkage
  f <- vapply(c(0, .25, .5, .75), function(s) get_compressed_angle(matrix(c(1, s, s, 1), 2)), 0)
  expect_true(all(diff(f) < 0))
})

test_that("invalid similarity matrices are rejected", {
  expect_error(get_compressed_angle(matrix(1:6, 2)), "square")
  expect_error(get_compressed_angle("a"), "square")
  expect_error(get_compressed_angle(matrix(c(1, NA, 0, 1), 2)), "missing")
})
