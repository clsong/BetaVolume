#' Geometric (hypervolume) beta diversity of a metacommunity
#'
#' @description
#' Computes the geometric beta diversity \eqn{\beta_{vol}} of Song et al.
#' (2025). The metacommunity matrix is embedded as a cloud of points (one point
#' per community or per species, whichever is more numerous) in a
#' \eqn{d}-dimensional space, \eqn{d = \min(\gamma, N)}{d = min(gamma, N)}, where \eqn{\gamma} is
#' the number of species and \eqn{N} the number of communities. Together with
#' the origin, the points span a convex polytope \eqn{P} and
#' \deqn{\beta_{vol} = d \times \mathrm{vol}(P)^{1/d}.}{beta_vol = d * vol(P)^(1/d).}
#' \eqn{\beta_{vol}} is an effective number of communities ranging from 0 (all
#' communities identical) to \eqn{d} (for presence/absence data: all
#' \eqn{2^d - 1} possible non-empty compositions are present). Unlike classic
#' measures, adding a distinct, species-rich community never decreases beta
#' diversity.
#'
#' @details
#' **Orientation.** The function is agnostic to whether communities are in rows
#' or columns: the matrix is transposed when necessary so that the smaller
#' dimension defines the axes (see [preprocess_meta_composition()]).
#'
#' **Why the origin?** An empty community (or a "ghost" species absent from
#' every community) must not change beta diversity (double-zero asymmetry), so
#' the origin is always part of the embedded object.
#'
#' **Estimation methods** (`hypervolume_method`):
#' \describe{
#'   \item{`"deterministic"` (alias `"exact"`)}{The exact volume of the convex
#'     hull, computed with Qhull through [geometry::convhulln()]. Computing
#'     convex hulls is exponential in \eqn{d} and is only practical up to
#'     about 15 dimensions.}
#'   \item{`"hyper_normal"`}{A fast approximation based on the ellipsoid
#'     spanned by the covariance matrix \eqn{\Sigma} of the embedded points
#'     (including the origin):
#'     \deqn{\beta \approx 4\, d\, \det(\Sigma)^{1/d}.}{beta ~ 4 d det(Sigma)^(1/d).}
#'     This is the variance-based geometric measure \eqn{\beta_{VAR}} (Box 1 of
#'     the paper) multiplied by \eqn{4 = 1/\mathrm{Var}(\mathrm{Bernoulli}(1/2))}{4 = 1/Var(Bernoulli(1/2))},
#'     which puts it on the same scale as \eqn{\beta_{vol}}: a
#'     presence/absence metacommunity containing every possible composition
#'     has covariance \eqn{\approx I/4}{~ I/4} and hence \eqn{\beta \approx d}{beta ~ d}. It
#'     scales to thousands of dimensions, accepts sparse matrices, and is
#'     the method used for large metacommunities in the paper. Note that,
#'     like a variance, it scales with the square of the data
#'     (\eqn{\beta(cZ) = c^2 \beta(Z)}) whereas \eqn{\beta_{vol}} scales
#'     linearly, so values from the two methods should not be mixed when
#'     comparing metacommunities with non-binary data.}
#' }
#' By default the exact method is used when \eqn{d \le}{d <=} `dim_threshold` (10)
#' and the hyper-normal approximation otherwise.
#'
#' **Species similarity and functional complementarity.** Species similarity
#' (Leinster & Cobbold 2012) acts as a linear transformation of the axes that
#' rescales the hypervolume by the normalised solid angle \eqn{\Omega(S)}
#' (Eq. 9 of the paper):
#' \eqn{\beta_{vol} = d\,(\Omega(S)\,\mathrm{vol}(P))^{1/d}}{beta_vol = d (Omega(S) vol(P))^(1/d)}. Pass
#' `compressed_angle = get_compressed_angle(S)` to apply it; see
#' [get_compressed_angle()].
#'
#' **Presence/absence duplicates.** Set `weights = TRUE` to use the weighted
#' embedding of duplicated compositions (Eq. 5 of the paper); see
#' [preprocess_meta_composition()].
#'
#' @inheritParams preprocess_meta_composition
#' @param preprocess Logical. If `TRUE` (default), `meta_composition` is passed
#'   through [preprocess_meta_composition()] (orientation, `weights`,
#'   `remove_unique`). Set to `FALSE` only if `meta_composition` is already a
#'   preprocessed embedding (rows = points, columns = axes).
#' @param hypervolume_method `"deterministic"` (exact convex hull; `"exact"`
#'   is accepted as an alias) or `"hyper_normal"` (fast ellipsoid
#'   approximation). If missing, chosen from `dim_threshold`.
#' @param compressed_angle Multiplicative factor \eqn{\Omega(S)^{1/d}}
#'   accounting for species similarity or functional complementarity, usually
#'   obtained from [get_compressed_angle()]. The default `1` corresponds to
#'   totally dissimilar species (identity similarity matrix).
#' @param dim_threshold When `hypervolume_method` is missing, the exact method
#'   is used if the embedding dimension is at most `dim_threshold` and the
#'   hyper-normal approximation otherwise. Default 10.
#' @param sparsity Logical. Is `meta_composition` a sparse matrix? Detected
#'   automatically; kept for backward compatibility.
#'
#' @return A single non-negative number: the geometric beta diversity.
#'
#' @references Song, C., Lu, M., Bennett, J. R., Gilbert, B., Fortin, M.-J., &
#'   Gonzalez, A. (2025). A geometric approach to beta diversity. *Ecological
#'   Monographs*, 95(1), e70008. \doi{10.1002/ecm.70008}
#'
#' Leinster, T., & Cobbold, C. A. (2012). Measuring diversity: the importance
#' of species similarity. *Ecology*, 93(3), 477-489. \doi{10.1890/10-2402.1}
#'
#' @seealso [preprocess_meta_composition()], [disentangle_individual_contribution()],
#'   [get_compressed_angle()], [beta_temporal_overlap()]
#' @export
#' @examples
#' # Fig. 1 of Song et al. (2025): two species A and B
#' meta_I <- rbind(c(1, 0), c(0, 1)) # {A}, {B}
#' meta_II <- rbind(c(1, 0), c(0, 1), c(1, 1)) # {A}, {B}, {A, B}
#' betavolume(meta_I) # 2 * sqrt(1/2) = 1.41
#' betavolume(meta_II) # 2 * sqrt(1) = 2, the maximum for two species
#'
#' # An empirical metacommunity (37 species x 10 plots)
#' betavolume(example_meta_composition)
#' betavolume(example_meta_composition, weights = TRUE)
#'
#' # Fast approximation for high-dimensional data
#' betavolume(example_meta_composition, hypervolume_method = "hyper_normal")
#'
#' # Sparse matrices are supported
#' sparse <- Matrix::Matrix(example_meta_composition, sparse = TRUE)
#' betavolume(sparse, hypervolume_method = "hyper_normal")
betavolume <- function(meta_composition,
                       preprocess = TRUE,
                       weights = FALSE,
                       remove_unique = FALSE,
                       hypervolume_method,
                       compressed_angle = 1,
                       dim_threshold = 10,
                       sparsity,
                       weight_scheme = c("paper", "legacy")) {
  check_flag(preprocess, "preprocess")
  if (!is.numeric(compressed_angle) || length(compressed_angle) != 1L ||
    !is.finite(compressed_angle) || compressed_angle < 0) {
    stop("`compressed_angle` must be a single non-negative number.", call. = FALSE)
  }
  weight_scheme <- match.arg(weight_scheme)

  if (preprocess) {
    points <- preprocess_meta_composition(meta_composition,
      weights = weights,
      remove_unique = remove_unique,
      weight_scheme = weight_scheme
    )
  } else {
    points <- as_composition_matrix(meta_composition)
  }
  if (!missing(sparsity)) {
    check_flag(sparsity, "sparsity")
    if (sparsity && !is_sparse(points)) {
      points <- Matrix::Matrix(points, sparse = TRUE)
    }
  }

  d <- ncol(points)
  if (nrow(points) == 0L) {
    return(0)
  }
  method <- resolve_method(
    if (missing(hypervolume_method)) NULL else hypervolume_method,
    d, dim_threshold
  )

  compressed_angle * beta_from_embedding(points, method)
}
