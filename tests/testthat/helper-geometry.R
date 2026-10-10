# Independent reference implementations used to validate BetaVolume.

# Area of the convex hull of 2-d points plus the origin, using base R's
# grDevices::chull() and the shoelace formula (no Qhull involved).
shoelace_hull_area <- function(points) {
  p <- rbind(as.matrix(points), 0)
  idx <- grDevices::chull(p[, 1], p[, 2])
  x <- p[idx, 1]
  y <- p[idx, 2]
  abs(sum(x * c(y[-1], y[1]) - c(x[-1], x[1]) * y)) / 2
}

# Volume of the convex hull of points plus the origin as the sum of the volumes
# of the simplices of a Delaunay triangulation (a different Qhull code path).
delaunay_volume <- function(points) {
  p <- unique(rbind(as.matrix(points), 0))
  d <- ncol(p)
  tri <- geometry::delaunayn(p, options = "Qt Qbb Qc Qz", output.options = FALSE)
  sum(apply(tri, 1, function(ix) {
    abs(det(t(p[ix[-1], , drop = FALSE]) - p[ix[1], ])) / factorial(d)
  }))
}

# Geometric beta diversity from a volume
beta_from_volume <- function(vol, d) d * vol^(1 / d)

# All 2^d presence/absence compositions (including the empty one)
all_compositions <- function(d) {
  as.matrix(expand.grid(rep(list(0:1), d)))
}

# Monte Carlo estimate of the fraction of R^d covered by the cone spanned by
# the columns of S: P(S^{-1} Z >= 0) for standard normal Z.
mc_cone_fraction <- function(S, n = 2e5) {
  z <- matrix(stats::rnorm(n * nrow(S)), nrow(S))
  coef <- solve(S, z)
  mean(colSums(coef >= 0) == nrow(S))
}

# Brute-force contributions (Eq. 7) by recomputing beta diversity without each
# unit (rows, or unique compositions when weighted).
brute_contribution <- function(m, weights = FALSE, method = "deterministic") {
  X <- preprocess_meta_composition(m, weights = weights)
  if (weights) {
    keys <- apply(X, 1, paste, collapse = "\r")
    units <- split(seq_len(nrow(X)), factor(keys, levels = unique(keys)))
  } else {
    units <- as.list(seq_len(nrow(X)))
  }
  b0 <- betavolume(X, preprocess = FALSE, hypervolume_method = method)
  loss <- vapply(units, function(u) {
    b0 - betavolume(X[-u, , drop = FALSE], preprocess = FALSE, hypervolume_method = method)
  }, numeric(1))
  unname(loss / sum(loss))
}
