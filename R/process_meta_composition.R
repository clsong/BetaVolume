#' Prepare a metacommunity matrix for the geometric embedding
#'
#' @description
#' Turns a metacommunity (community composition) matrix into the point cloud
#' whose hypervolume defines geometric beta diversity (Song et al. 2025):
#'
#' 1. optionally removes duplicated rows and duplicated columns
#'    (`remove_unique = TRUE`);
#' 2. orients the matrix so that **rows are the embedded points and columns the
#'    axes** of the embedding space. The smaller of the two dimensions
#'    (`min(number of species, number of communities)`) is the embedding
#'    dimension \eqn{d}, so the matrix is transposed when it has fewer rows than
#'    columns;
#' 3. optionally applies the weighted embedding for duplicated compositions
#'    (`weights = TRUE`, Song et al. 2025, Eq. 5 and Fig. 4).
#'
#' You rarely need to call this function directly: [betavolume()],
#' [disentangle_individual_contribution()] and [beta_temporal_overlap()] call it
#' for you unless `preprocess = FALSE`.
#'
#' @section Weighted embedding of duplicated compositions:
#' With presence/absence data several communities often share exactly the same
#' composition and are therefore embedded at the same point. With
#' `weights = TRUE` all-zero rows are dropped and each row is multiplied by a
#' weight proportional to the number of times its composition occurs. With
#' `weight_scheme = "paper"` (the default) the weight of a composition that
#' occurs \eqn{n_k} times is
#' \deqn{w_k = \frac{m \, n_k}{\sum_{l=1}^{m} n_l},}{w_k = m n_k / sum(n_l),}
#' where \eqn{m} is the number of unique compositions (Eq. 5 of the paper). The
#' weights average to one over the unique compositions, so a metacommunity
#' without duplicates is left unchanged. For the six-community example of
#' Fig. 4 of the paper the weights are \eqn{3 \times 1/6}{3 x 1/6}, \eqn{3 \times 2/6}{3 x 2/6}
#' and \eqn{3 \times 3/6}{3 x 3/6}.
#'
#' `weight_scheme = "legacy"` reproduces BetaVolume <= 0.1.0, which used
#' \eqn{w_k = N n_k / \sum_l n_l^2}{w_k = N n_k / sum(n_l^2)} with \eqn{N} the number of non-empty rows.
#' Both schemes are proportional to \eqn{n_k} and therefore give **identical
#' relative contributions** ([disentangle_individual_contribution()]), but the
#' legacy weights are smaller by the factor \eqn{m \sum_l n_l^2 / N^2 \ge 1}{m sum(n_l^2) / N^2 >= 1},
#' so `betavolume(weights = TRUE)` values differ. Use `"legacy"` only to
#' reproduce results computed with older versions of the package.
#'
#' @param meta_composition A numeric matrix (or data frame, or sparse
#'   [Matrix::Matrix]) with communities in rows and species in columns (or the
#'   reverse -- orientation does not matter). Entries are presence/absence
#'   (0/1) or any other non-negative measure of species importance such as
#'   abundance or biomass.
#' @param weights Logical. Apply the weighted embedding of duplicated
#'   compositions (Eq. 5 of Song et al. 2025; see the section "Weighted
#'   embedding of duplicated compositions" in [preprocess_meta_composition()])?
#'   Default `FALSE`.
#' @param remove_unique Logical. Remove duplicated rows and duplicated columns
#'   before embedding? Default `FALSE`. (Despite the historical name, this
#'   keeps one copy of each *unique* row/column.)
#' @param weight_scheme Either `"paper"` (default, Eq. 5 of Song et al. 2025) or
#'   `"legacy"` (weights used by BetaVolume <= 0.1.0). Only used when
#'   `weights = TRUE`.
#'
#' @return A numeric matrix (sparse if the input was sparse) whose rows are the
#'   embedded points and whose columns are the axes, i.e. with
#'   `nrow >= ncol`. Row names are kept; if the input had none, rows are named
#'   by their index in the (possibly transposed) input so that they can be
#'   matched to the original data after empty rows have been dropped.
#'
#' @references Song, C., Lu, M., Bennett, J. R., Gilbert, B., Fortin, M.-J., &
#'   Gonzalez, A. (2025). A geometric approach to beta diversity. *Ecological
#'   Monographs*, 95(1), e70008. \doi{10.1002/ecm.70008}
#'
#' @seealso [betavolume()]
#' @export
#' @examples
#' # Fig. 4 of Song et al. (2025): two species, six communities
#' fig4 <- rbind(
#'   c(1, 0),
#'   c(0, 1), c(0, 1),
#'   c(1, 1), c(1, 1), c(1, 1)
#' )
#' colnames(fig4) <- c("A", "B")
#' preprocess_meta_composition(fig4, weights = TRUE)
#'
#' # Wide matrices are transposed so that rows are the embedded points
#' dim(preprocess_meta_composition(t(fig4)))
preprocess_meta_composition <- function(meta_composition,
                                        weights = FALSE,
                                        remove_unique = FALSE,
                                        weight_scheme = c("paper", "legacy")) {
  check_flag(weights, "weights")
  check_flag(remove_unique, "remove_unique")
  weight_scheme <- match.arg(weight_scheme)
  x <- as_composition_matrix(meta_composition)

  if (remove_unique) {
    x <- x[!duplicated(row_keys(x)), , drop = FALSE]
    x <- x[, !duplicated(row_keys(transpose(x))), drop = FALSE]
  }

  if (nrow(x) < ncol(x)) {
    x <- transpose(x)
  }

  if (is.null(rownames(x))) {
    rownames(x) <- as.character(seq_len(nrow(x)))
  }

  if (weights) {
    x <- weight_composition(x, weight_scheme)
  }
  x
}

#' Weighted embedding of duplicated compositions (Song et al. 2025, Eq. 5)
#' @noRd
weight_composition <- function(x, weight_scheme = "paper") {
  nonzero <- if (is_sparse(x)) Matrix::rowSums(x != 0) > 0 else rowSums(x != 0) > 0
  x <- x[nonzero, , drop = FALSE]
  if (nrow(x) == 0L) {
    return(x)
  }
  groups <- row_groups(x)
  n_k <- groups$counts[groups$group] # multiplicity of each row's composition
  n_total <- nrow(x)
  m_unique <- length(groups$counts)
  w <- switch(weight_scheme,
    paper = m_unique * n_k / n_total,
    legacy = n_total * n_k / sum(groups$counts^2)
  )
  if (is_sparse(x)) {
    rn <- rownames(x)
    x <- Matrix::Diagonal(x = w) %*% x
    rownames(x) <- rn
    x
  } else {
    x * w
  }
}
