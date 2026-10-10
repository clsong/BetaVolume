# Fig. 5 of Song et al. (2025): rows = communities, columns = species A and B
t0 <- rbind(c(0, 1), c(1, 1)) # community 1 = {B}, community 2 = {A, B}
t1 <- rbind(c(0, 1), c(1, 0)) # community 2 loses species B
t2 <- rbind(c(1, 0), c(0, 1)) # the two communities swap compositions

test_that("Fig. 5: asynchronous change t -> t+1", {
  # vol(P_t) = 1/2, vol(P_t+1 n P_t) = 1/4  =>  (1/4 / 1/2)^(1/2) = 0.71
  expect_equal(beta_temporal_overlap(t0, t1), sqrt(1 / 2))
  expect_equal(beta_temporal_overlap(t0, t1, oriented = TRUE), sqrt(1 / 2))
})

test_that("Fig. 5: synchronous change t+1 -> t+2 (opposite orientation)", {
  # identical shapes: unoriented overlap 1; oriented: (0.5 - (-0.5)) / 0.5 => sqrt(2)
  expect_equal(beta_temporal_overlap(t1, t2), 1)
  expect_equal(beta_temporal_overlap(t1, t2, oriented = TRUE), sqrt(2))
})

test_that("identical and nested metacommunities", {
  set.seed(31)
  m <- matrix(stats::runif(20 * 3), 20, 3)
  expect_equal(beta_temporal_overlap(m, m), 1)
  # P2 = P1 / 2 is contained in P1 and has 2^-d of its volume
  expect_equal(beta_temporal_overlap(m, m / 2), 0.5, tolerance = 1e-8)
  # P2 = 2 P1 contains P1
  expect_equal(beta_temporal_overlap(m, 2 * m), 1, tolerance = 1e-8)
})

test_that("the overlap of hulls sharing only a face is 0", {
  a <- rbind(c(1, 0), c(1, 1)) # triangle below the diagonal
  b <- rbind(c(0, 1), c(1, 1)) # triangle above the diagonal
  expect_equal(beta_temporal_overlap(a, b), 0)
})

test_that("2-d overlap matches an analytic intersection", {
  # Unit square vs the triangle (0,0), (1,0), (0,1): overlap 1/2 of the square
  sq <- rbind(c(1, 0), c(0, 1), c(1, 1))
  tri <- rbind(c(1, 0), c(0, 1), c(0, 0))
  expect_equal(beta_temporal_overlap(sq, tri), sqrt(1 / 2), tolerance = 1e-10)
  expect_equal(beta_temporal_overlap(tri, sq), 1, tolerance = 1e-10)
})

test_that("invalid inputs are rejected", {
  expect_error(beta_temporal_overlap(t0, cbind(t0, 1)), "same dimensions")
  expect_error(beta_temporal_overlap(t0, t1, remove_unique = TRUE), "remove_unique")
  m3 <- rbind(c(1, 0), c(0, 1), c(1, 1))
  expect_error(beta_temporal_overlap(m3, m3, oriented = TRUE), "simplex")
  expect_warning(res <- beta_temporal_overlap(rbind(c(1, 1), c(2, 2)), t1), "zero volume")
  expect_true(is.na(res))
})
