#' Contribution of each community (or species) to geometric beta diversity
#'
#' @description
#' Quantifies how much each embedded point -- a community or a species,
#' whichever is the more numerous dimension of the metacommunity matrix --
#' contributes to the geometric beta diversity, as the normalised loss of
#' beta diversity when that point is removed (Song et al. 2025, Eq. 7 and
#' Fig. 6):
#' \deqn{C_i = \frac{\mathrm{vol}(P_0)^{1/d} - \mathrm{vol}(P_{\setminus i})^{1/d}}
#'   {\sum_j \left(\mathrm{vol}(P_0)^{1/d} - \mathrm{vol}(P_{\setminus j})^{1/d}\right)},}{
#'   C_i = (vol(P_0)^(1/d) - vol(P_-i)^(1/d)) / sum_j (vol(P_0)^(1/d) - vol(P_-j)^(1/d)),}
#' where \eqn{P_0} is the geometric object of the full metacommunity and
#' \eqn{P_{\setminus i}}{P_-i} the object without point \eqn{i}.
#'
#' @details
#' **Unit of removal.** Without weighting (`weights = FALSE`) each row of the
#' embedding is removed in turn. If a composition occurs several times,
#' removing one copy leaves the geometric object unchanged, so each copy
#' contributes 0 (only communities with a unique composition can contribute).
#' With the weighted embedding (`weights = TRUE`), duplicated compositions are
#' merged into a single weighted point (Eq. 5), and each *unique composition*
#' is removed in turn, as in Fig. 6 of the paper. The result then has one
#' value per unique composition, named by the rows that share it (for example
#' `"2, 3"`). The weights are computed once on the full metacommunity and are
#' not recomputed after a composition is removed. Both weighting schemes are
#' proportional to the multiplicity of each composition and \eqn{C_i} is scale
#' invariant, so `"paper"` and `"legacy"` weights give identical contributions.
#'
#' **Sign.** With the exact convex-hull method every contribution is
#' non-negative, because removing a point never increases the volume of a
#' convex hull. Only vertices of the hull contribute; points inside the hull
#' contribute 0. With the `"hyper_normal"` approximation, removing a point close
#' to the centroid can *increase* the determinant of the covariance matrix, so
#' small negative contributions are possible.
#'
#' If the embedding is degenerate (zero beta diversity) the contributions are
#' undefined and `NA` is returned with a warning. If beta diversity is positive
#' but no single point is essential (every vertex of the hull occurs more than
#' once and `weights = FALSE`), all contributions are 0 and a warning suggests
#' `weights = TRUE`.
#'
#' @section Performance:
#' A naive implementation recomputes one convex hull per point. Here, for the
#' exact method only the vertices of the hull that occur once are recomputed.
#' For the hyper-normal method all leave-one-out determinants are obtained at
#' once from a low-rank downdate of the covariance matrix (matrix determinant
#' lemma) at the cost of a single Cholesky factorisation, which is two to
#' three orders of magnitude faster than recomputing every determinant. The
#' exact hulls can also be recomputed in parallel with `cores > 1`.
#'
#' @inheritParams betavolume
#' @param cores Number of cores used to recompute the convex hulls of the
#'   exact method in parallel (forked processes via [parallel::mclapply()];
#'   ignored on Windows). Defaults to `getOption("BetaVolume.cores", 1L)`.
#'
#' @return A named numeric vector of contributions summing to one, with one
#'   value per embedded point (row of
#'   `preprocess_meta_composition(meta_composition, ...)`), or one value per
#'   unique composition when `weights = TRUE`. Names are the row names of the
#'   embedding (community or species names, or their indices).
#'
#' @references Song, C., Lu, M., Bennett, J. R., Gilbert, B., Fortin, M.-J., &
#'   Gonzalez, A. (2025). A geometric approach to beta diversity. *Ecological
#'   Monographs*, 95(1), e70008. \doi{10.1002/ecm.70008}
#'
#' @seealso [betavolume()]
#' @export
#' @examples
#' # Fig. 4 and 6 of Song et al. (2025): six communities, three unique
#' # compositions {A}, {B} (twice) and {A, B} (three times)
#' fig4 <- rbind(c(1, 0), c(0, 1), c(0, 1), c(1, 1), c(1, 1), c(1, 1))
#' round(disentangle_individual_contribution(fig4, weights = TRUE), 2)
#'
#' # The same result from the weighted matrix shown in Fig. 6A
#' fig6 <- rbind(c(0.5, 0), c(0, 1), c(1.5, 1.5))
#' round(disentangle_individual_contribution(fig6), 2)
#'
#' # Contributions of the 37 species of the example metacommunity
#' contrib <- disentangle_individual_contribution(example_meta_composition)
#' head(sort(contrib, decreasing = TRUE))
disentangle_individual_contribution <- function(meta_composition,
                                                weights = FALSE,
                                                remove_unique = FALSE,
                                                hypervolume_method,
                                                dim_threshold = 10,
                                                preprocess = TRUE,
                                                weight_scheme = c("paper", "legacy"),
                                                cores = getOption("BetaVolume.cores", 1L)) {
  check_flag(preprocess, "preprocess")
  if (!is.numeric(cores) || length(cores) != 1L || is.na(cores) || cores < 1) {
    stop("`cores` must be a positive integer.", call. = FALSE)
  }
  check_flag(weights, "weights")
  weight_scheme <- match.arg(weight_scheme)
  if (preprocess) {
    points <- preprocess_meta_composition(meta_composition,
      weights = weights,
      remove_unique = remove_unique,
      weight_scheme = weight_scheme
    )
  } else {
    points <- as_composition_matrix(meta_composition)
    if (is.null(rownames(points))) rownames(points) <- as.character(seq_len(nrow(points)))
  }
  ids <- rownames(points)
  d <- ncol(points)
  method <- resolve_method(
    if (missing(hypervolume_method)) NULL else hypervolume_method,
    d, dim_threshold
  )

  if (nrow(points) == 0L) {
    return(stats::setNames(numeric(0), character(0)))
  }

  # Units removed one at a time: rows, or unique compositions when weighted.
  if (weights && preprocess) {
    groups <- row_groups(points)
    group <- groups$group
    labels <- vapply(split(ids, factor(group, levels = seq_along(groups$first))),
      paste, character(1),
      collapse = ", "
    )
  } else {
    group <- seq_len(nrow(points))
    labels <- ids
  }

  loss <- switch(method,
    deterministic = {
      if (weights && preprocess) {
        exact_loo_loss(as.matrix(points[groups$first, , drop = FALSE]), cores)
      } else {
        exact_loo_loss(as.matrix(points), cores)
      }
    },
    hyper_normal = hyper_normal_loo_loss(points, group)
  )

  labels <- unname(labels)
  if (isTRUE(attr(loss, "degenerate"))) {
    warning("The metacommunity has zero beta diversity; contributions are undefined.",
      call. = FALSE
    )
    return(stats::setNames(rep(NA_real_, length(loss)), labels))
  }
  total <- sum(loss)
  if (!is.finite(total) || total <= 1e-10 * attr(loss, "scale")) {
    warning("Removing any single point leaves beta diversity unchanged (every vertex of ",
      "the hull occurs more than once), so all contributions are 0. Use ",
      "`weights = TRUE` to treat duplicated compositions as single units.",
      call. = FALSE
    )
    return(stats::setNames(numeric(length(loss)), labels))
  }
  stats::setNames(as.vector(loss) / total, labels)
}

#' Leave-one-out loss of vol^(1/d) for the exact convex hull.
#'
#' Only hull vertices whose composition occurs once can change the hull when
#' removed, so only those require a new hull computation.
#' @noRd
exact_loo_loss <- function(points, cores = 1L) {
  n <- nrow(points)
  d <- ncol(points)
  loss <- numeric(n)

  keys <- row_keys(points)
  is_origin <- rowSums(points != 0) == 0
  # Unique non-origin points (the origin is always part of P).
  keep <- !duplicated(keys) & !is_origin
  multiplicity <- tabulate(match(keys, keys[keep]), nbins = sum(keep))
  uniq <- points[keep, , drop = FALSE]
  owner <- which(keep) # row of `points` represented by each unique point

  all_pts <- rbind(uniq, 0) # origin is the last row
  hull <- qhull(all_pts)
  if (is.null(hull)) {
    attr(loss, "degenerate") <- TRUE
    return(loss)
  }
  v0 <- hull$vol^(1 / d)
  attr(loss, "scale") <- v0
  vertices <- sort(unique(as.vector(hull$hull)))
  vertices <- vertices[vertices <= nrow(uniq)] # drop the origin
  candidates <- vertices[multiplicity[vertices] == 1L]

  loo_root_volume <- function(k) {
    h <- qhull(all_pts[-k, , drop = FALSE])
    if (is.null(h)) 0 else h$vol^(1 / d)
  }
  vk <- if (cores > 1L && length(candidates) > 1L && .Platform$OS.type == "unix") {
    res <- parallel::mclapply(candidates, loo_root_volume, mc.cores = as.integer(cores))
    failed <- vapply(res, inherits, logical(1), what = "try-error")
    if (any(failed)) {
      stop("Parallel hull computation failed: ", as.character(res[[which(failed)[1]]]), call. = FALSE)
    }
    unlist(res)
  } else {
    vapply(candidates, loo_root_volume, numeric(1))
  }
  loss[owner[candidates]] <- pmax(v0 - vk, 0)
  loss
}

#' Leave-one-group-out loss for the hyper-normal approximation.
#'
#' Let X contain the points plus the origin (n rows), with mean mu, covariance
#' C and scatter matrix M = (n - 1) C. Removing a group of n_k identical points
#' x is a rank-one downdate of the scatter matrix,
#'   M_{-k} = M - n n_k / (n - n_k) * dd',   d = x - mu,
#' so by the matrix determinant lemma
#'   det(M_{-k}) = det(M) * (1 - n n_k / ((n - n_k)(n - 1)) * d' C^{-1} d),
#' and all leave-one-out determinants cost a single Cholesky factorisation.
#' With singleton groups this is the usual leave-one-out downdate.
#'
#' @param group Integer vector giving, for each row of `points`, the unit it
#'   belongs to (rows of a unit must be identical).
#' @return The loss of det(C)^(1/d) for each unit (the factor 4 d of the
#'   hyper-normal estimator cancels in the normalised contributions).
#' @noRd
hyper_normal_loo_loss <- function(points, group = seq_len(nrow(points))) {
  d <- ncol(points)
  n <- nrow(points) + 1L # + origin
  n_units <- max(group)
  first <- match(seq_len(n_units), group)
  n_k <- tabulate(group, nbins = n_units)

  cov <- embedding_covariance(points)
  ld <- logdet_psd(cov)
  r <- suppressWarnings(tryCatch(chol(cov), error = function(e) NULL))
  if (is.null(r) || !is.finite(ld)) {
    return(structure(rep(0, n_units), degenerate = TRUE))
  }
  mu <- if (is_sparse(points)) Matrix::colSums(points) / n else colSums(points) / n
  x <- as.matrix(points[first, , drop = FALSE])
  delta <- x - rep(mu, each = n_units)
  # Mahalanobis distances d' C^{-1} d via a triangular solve: C = R'R
  z <- backsolve(r, t(delta), transpose = TRUE)
  maha <- colSums(z^2)

  beta0 <- exp(ld / d)
  remaining <- n - n_k # points left, origin included
  factor <- 1 - n * n_k / ((n - n_k) * (n - 1)) * maha
  logdet_scatter <- ld + d * log(n - 1)
  ok <- remaining >= 2L & factor > sqrt(.Machine$double.eps)
  beta_loo <- numeric(n_units)
  beta_loo[ok] <- exp((logdet_scatter + log(factor[ok]) - d * log(remaining[ok] - 1)) / d)
  structure(beta0 - beta_loo, scale = beta0)
}
