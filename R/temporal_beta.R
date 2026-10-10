#' Temporal change of beta diversity as hypervolume overlap
#'
#' @description
#' Measures how the geometric embedding of a metacommunity changes between two
#' times as the (oriented) overlap of the two hypervolumes (Song et al. 2025,
#' Eq. 6 and Fig. 5):
#' \deqn{\beta_{t \to t+1} = \left(\frac{\mathrm{vol}(P_{t+1} \cap P_t)}{\mathrm{vol}(P_t)}\right)^{1/d}.}{
#'   beta_{t -> t+1} = (vol(P_{t+1} intersect P_t) / vol(P_t))^(1/d).}
#'
#' @details
#' **Unoriented overlap** (`oriented = FALSE`, default). The ratio of the
#' volume of the intersection of the two convex hulls (both including the
#' origin) to the volume of the hull at time 1. It lies in \eqn{[0, 1]}: 1 if
#' the metacommunity at time 2 covers the one at time 1, 0 if the two
#' embeddings only share a lower-dimensional face.
#'
#' **Oriented overlap** (`oriented = TRUE`). The paper assigns an orientation
#' to the embedding -- from the origin to community 1, to community 2, and so
#' on -- to distinguish synchronous from asynchronous changes. The orientation
#' (sign of the volume) is defined when each metacommunity is a simplex, i.e.
#' when it has exactly \eqn{d} embedded points (for example \eqn{N} communities
#' and \eqn{N} species, as in Fig. 5), given in the same order at both times. If
#' both embeddings have the same orientation, the overlap is the volume of the
#' intersection, as above; if the orientations are opposite, the overlap is the
#' signed difference \eqn{|\mathrm{vol}(P_t)| + |\mathrm{vol}(P_{t+1})|}{|vol(P_t)| + |vol(P_{t+1})|}.
#' Values below 1 indicate asynchronous changes (or synchronous changes in the
#' same direction), values above 1 synchronous changes in opposite directions.
#'
#' Both metacommunities must describe the same species and communities, i.e.
#' have the same dimensions with rows and columns in the same order.
#'
#' @param meta_composition_1,meta_composition_2 Metacommunity matrices at time
#'   1 and time 2 (see [betavolume()]), with identical dimensions.
#' @param preprocess Logical. Pass both matrices through
#'   [preprocess_meta_composition()]? Default `TRUE`.
#' @param oriented Logical. Use the oriented overlap (see Details)? Default
#'   `FALSE`.
#' @param ... Further arguments passed to [preprocess_meta_composition()]
#'   (`weights`, `weight_scheme`). `remove_unique = TRUE` is not allowed
#'   because it could remove different species at the two times and misalign
#'   the axes.
#'
#' @return A single number. `NA` (with a warning) if the metacommunity at time
#'   1 has zero beta diversity, in which case the ratio is undefined.
#'
#' @references Song, C., Lu, M., Bennett, J. R., Gilbert, B., Fortin, M.-J., &
#'   Gonzalez, A. (2025). A geometric approach to beta diversity. *Ecological
#'   Monographs*, 95(1), e70008. \doi{10.1002/ecm.70008}
#'
#' @seealso [betavolume()]
#' @export
#' @examples
#' # Fig. 5 of Song et al. (2025); rows = communities, columns = species A, B
#' t0 <- rbind(c(0, 1), c(1, 1)) # community 1 = {B}, community 2 = {A, B}
#' t1 <- rbind(c(0, 1), c(1, 0)) # community 2 loses species B
#' t2 <- rbind(c(1, 0), c(0, 1)) # communities 1 and 2 swap compositions
#'
#' beta_temporal_overlap(t0, t1) # (1/4 / 1/2)^(1/2) = 0.71
#' beta_temporal_overlap(t1, t2) # identical shapes: 1
#' beta_temporal_overlap(t1, t2, oriented = TRUE) # opposite orientation: 1.41
beta_temporal_overlap <- function(meta_composition_1,
                                  meta_composition_2,
                                  preprocess = TRUE,
                                  oriented = FALSE,
                                  ...) {
  check_flag(preprocess, "preprocess")
  check_flag(oriented, "oriented")
  m1 <- as_composition_matrix(meta_composition_1, "meta_composition_1")
  m2 <- as_composition_matrix(meta_composition_2, "meta_composition_2")
  if (!identical(dim(m1), dim(m2))) {
    stop("`meta_composition_1` and `meta_composition_2` must have the same dimensions ",
      "(the same species and communities in the same order).",
      call. = FALSE
    )
  }
  if (preprocess) {
    dots <- list(...)
    if (isTRUE(dots$remove_unique)) {
      stop("`remove_unique = TRUE` is not supported by beta_temporal_overlap() ",
        "because it can misalign the axes of the two embeddings.",
        call. = FALSE
      )
    }
    m1 <- preprocess_meta_composition(m1, ...)
    m2 <- preprocess_meta_composition(m2, ...)
  }
  x1 <- as.matrix(m1)
  x2 <- as.matrix(m2)
  if (ncol(x1) != ncol(x2)) {
    stop("The two embeddings have different dimensions.", call. = FALSE)
  }
  d <- ncol(x1)
  p1 <- hull_points(x1)
  p2 <- hull_points(x2)

  if (!oriented) {
    v1 <- hull_volume(x1)
    if (v1 <= 0) {
      warning("The metacommunity at time 1 has zero volume; the overlap is undefined.",
        call. = FALSE
      )
      return(NA_real_)
    }
    v12 <- intersection_volume(p1, p2)
    return(min(v12 / v1, 1)^(1 / d))
  }

  if (nrow(x1) != d || nrow(x2) != d) {
    stop("The oriented overlap requires each embedding to be a simplex: exactly ", d,
      " embedded points (non-empty communities) in ", d, " dimensions.",
      call. = FALSE
    )
  }
  s1 <- det(x1) / factorial(d)
  s2 <- det(x2) / factorial(d)
  if (s1 == 0) {
    warning("The metacommunity at time 1 has zero volume; the overlap is undefined.",
      call. = FALSE
    )
    return(NA_real_)
  }
  overlap <- if (s2 == 0 || sign(s1) == sign(s2)) {
    intersection_volume(p1, p2)
  } else {
    abs(s1) + abs(s2)
  }
  (overlap / abs(s1))^(1 / d)
}

#' Volume of the intersection of two convex hulls given by their points.
#' @noRd
intersection_volume <- function(p1, p2) {
  d <- ncol(p1)
  if (d == 1L) {
    lo <- max(min(p1), min(p2))
    hi <- min(max(p1), max(p2))
    return(max(hi - lo, 0))
  }
  if (!is_full_dimensional(p1) || !is_full_dimensional(p2)) {
    return(0)
  }
  res <- tryCatch(
    geometry::intersectn(p1, p2, options = "Qt", return.chs = FALSE),
    error = function(e) e
  )
  if (inherits(res, "error")) {
    if (grepl("QH6154|QH6013|QH6023|flat|feasible", conditionMessage(res))) {
      return(0) # intersection is lower dimensional
    }
    stop("Failed to intersect the two hulls: ", conditionMessage(res), call. = FALSE)
  }
  res$ch$vol
}
