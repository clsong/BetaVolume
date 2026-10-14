# Hypervolume estimators ------------------------------------------------------
#
# Internal workhorses. Every function takes the *preprocessed* embedding: a
# matrix whose rows are the embedded points and whose d columns are the axes.
# The origin is added here (it is always part of the geometric object P, see
# Song et al. 2025, Eq. 3), so callers must not add it themselves.

#' Qhull options used for every exact convex hull.
#'
#' `"Qt"` triangulates the output. We deliberately do *not* use `"Tv"`
#' (verification), the default of `geometry::convhulln()`: on degenerate
#' presence/absence data with many coplanar points `"Tv"` can abort with an
#' internal "qh_checkpolygon" error even though the hull itself is fine. The
#' legacy implementation caught such errors and silently returned 0.
#'
#' `"Q0"` turns off Qhull's pre-merging of nearly coplanar facets (in 5-d and
#' higher Qhull otherwise adds `"Qx"`, exact pre-merges). Presence/absence
#' points have many coplanar facets, and pre-merging them dominates the run
#' time; without it hulls are about four times faster and their volumes agree
#' with the default options to about 1e-13. If Qhull reports an error with
#' `"Q0"`, the hull is recomputed with the default options
#' (`qhull_fallback_options`).
#' @noRd
qhull_options <- "Qt Q0"
qhull_fallback_options <- "Qt"

#' Prepare a point cloud for Qhull: densify, add the origin, drop duplicates.
#' @noRd
hull_points <- function(points) {
  points <- as.matrix(points)
  points <- rbind(points, 0)
  points[!duplicated(row_keys(points)), , drop = FALSE]
}

#' Is the convex hull of `points` (which contain the origin) full dimensional?
#'
#' Because the origin is one of the points, the affine hull is full
#' dimensional iff the points span R^d linearly, i.e. iff `rank(points) == d`.
#' @noRd
is_full_dimensional <- function(points, tol = 1e-10) {
  d <- ncol(points)
  if (nrow(points) < d + 1L) {
    return(FALSE)
  }
  scale <- max(abs(points))
  if (scale == 0) {
    return(FALSE)
  }
  qr(points / scale, tol = tol)$rank >= d
}

#' Run Qhull and return the hull (with volume), handling degenerate input.
#'
#' @return `NULL` when the hull is degenerate (zero volume), otherwise the
#'   output of `geometry::convhulln()` with `$hull` and `$vol`.
#' @noRd
qhull <- function(points) {
  if (ncol(points) == 1L) {
    # Qhull needs at least two dimensions; in 1-d the hull is an interval.
    rng <- range(points)
    if (diff(rng) <= 0) {
      return(NULL)
    }
    lo <- which.min(points[, 1L])
    hi <- which.max(points[, 1L])
    return(list(hull = matrix(c(lo, hi), ncol = 1L), vol = diff(rng)))
  }
  if (!is_full_dimensional(points)) {
    return(NULL)
  }
  run <- function(options) {
    geometry::convhulln(points, options = options, output.options = "FA")
  }
  on_error <- function(e) {
    msg <- conditionMessage(e)
    if (grepl("QH6154|QH6013|flat|less than|not enough points", msg)) {
      # Precision-level flat input: the volume is numerically zero.
      return(NULL)
    }
    # Retry once with joggled input before giving up loudly.
    tryCatch(run(paste(qhull_fallback_options, "QJ")), error = function(e2) {
      stop("Qhull failed to compute the convex hull: ", msg, call. = FALSE)
    })
  }
  tryCatch(run(qhull_options), error = function(e) {
    # Without pre-merging Qhull can stop on a precision error that the
    # default options handle: retry with those first.
    tryCatch(run(qhull_fallback_options), error = on_error)
  })
}

#' Exact volume of the convex hull of `points` together with the origin.
#' @noRd
hull_volume <- function(points) {
  h <- qhull(hull_points(points))
  if (is.null(h)) 0 else h$vol
}

#' Covariance matrix of the embedded points plus the origin.
#'
#' Works for dense and sparse (`Matrix`) inputs. For sparse inputs the
#' covariance is formed as `(X'X - n * m m') / (n - 1)` so that the data are
#' never densified (only the d x d result is dense). Equals
#' `stats::cov(rbind(points, 0))`.
#' @noRd
embedding_covariance <- function(points) {
  n <- nrow(points) + 1L # + origin
  if (is_sparse(points)) {
    mu <- Matrix::colSums(points) / n
    xtx <- as.matrix(Matrix::crossprod(points))
    cov <- (xtx - n * tcrossprod(mu)) / (n - 1)
  } else {
    # Centred cross-product (BLAS) with the origin row folded in analytically:
    # the origin contributes (0 - mu)(0 - mu)' = mu mu'.
    mu <- colSums(points) / n
    centred <- points - rep(mu, each = nrow(points))
    cov <- (crossprod(centred) + tcrossprod(mu)) / (n - 1)
  }
  dimnames(cov) <- NULL
  (cov + t(cov)) / 2
}

#' log-determinant of a symmetric positive semi-definite matrix.
#'
#' Returns `-Inf` when the matrix is singular (not numerically positive
#' definite). A Cholesky factorisation is about an order of magnitude faster
#' than the eigen-decomposition used previously and avoids `NaN`s produced by
#' tiny negative eigenvalues.
#' @noRd
logdet_psd <- function(m) {
  if (length(m) == 0L) {
    return(0)
  }
  r <- suppressWarnings(tryCatch(chol(m), error = function(e) NULL))
  if (is.null(r)) {
    return(-Inf)
  }
  dg <- diag(r)
  if (any(dg <= sqrt(.Machine$double.eps) * max(dg))) {
    return(-Inf)
  }
  2 * sum(log(dg))
}

#' Calibration constant of the hyper-normal approximation.
#'
#' `1 / Var(Bernoulli(1/2)) = 4`: with it, a presence/absence metacommunity that
#' contains every one of the 2^d possible compositions (covariance I/4) has a
#' hyper-normal beta diversity of (approximately) d, the maximum of beta_vol.
#' @noRd
hyper_normal_constant <- 4

#' Hyper-normal (ellipsoid) beta diversity: 4 * d * det(Cov)^(1/d)
#' @noRd
hyper_normal_beta <- function(points) {
  d <- ncol(points)
  ld <- logdet_psd(embedding_covariance(points))
  if (!is.finite(ld)) {
    return(0)
  }
  hyper_normal_constant * d * exp(ld / d)
}

#' Exact beta diversity d * vol(P)^(1/d)
#' @noRd
exact_beta <- function(points) {
  d <- ncol(points)
  d * hull_volume(points)^(1 / d)
}

#' Dispatch on the method
#' @noRd
beta_from_embedding <- function(points, method) {
  switch(method,
    deterministic = exact_beta(points),
    hyper_normal = hyper_normal_beta(points)
  )
}
