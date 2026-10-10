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
#' weight proportional to the number of times \eqn{n_k} its composition
#' occurs. Two normalisations are available, with \eqn{N} the number of
#' non-empty communities and \eqn{m} the number of unique compositions:
#' \describe{
#'   \item{`"legacy"` (default)}{\eqn{w_k = N n_k / \sum_l n_l^2}{w_k = N n_k / sum(n_l^2)}.
#'     The weights average to one over communities. This is the normalisation
#'     used for all analyses of Song et al. (2025), including Figs. 9 and 10,
#'     and in practice it keeps \eqn{\beta_{vol}} within its range
#'     \eqn{[0, d]}.}
#'   \item{`"paper"`}{\eqn{w_k = m n_k / N}{w_k = m n_k / N}, as printed in
#'     Eq. 5 and Fig. 4 of the paper. The weights average to one over unique
#'     compositions. Because they can exceed one substantially,
#'     \eqn{\beta_{vol}} can exceed its nominal maximum \eqn{d}: the Fig. 4
#'     example gives 2.12 for two species, and about one in ten random
#'     presence/absence metacommunities exceeds \eqn{d}.}
#' }
#' Both leave a metacommunity without duplicates unchanged. They are
#' proportional to each other, differing by the factor
#' \eqn{m \sum_l n_l^2 / N^2 \ge 1}{m sum(n_l^2) / N^2 >= 1}, so they give
#' **identical relative contributions**
#' ([disentangle_individual_contribution()]) but different
#' `betavolume(weights = TRUE)` values.
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
#' @param weight_scheme Normalisation of the weights: `"legacy"` (default; the
#'   weights used for the analyses of Song et al. 2025) or `"paper"` (Eq. 5 as
#'   printed in the paper). See the section on the weighted embedding. Only
#'   used when `weights = TRUE`.
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
#' preprocess_meta_composition(fig4, weights = TRUE) # default "legacy" weights
#' preprocess_meta_composition(fig4, weights = TRUE, weight_scheme = "paper") # Eq. 5
#'
#' # Wide matrices are transposed so that rows are the embedded points
#' dim(preprocess_meta_composition(t(fig4)))
preprocess_meta_composition <- function(meta_composition,
                                        weights = FALSE,
                                        remove_unique = FALSE,
                                        weight_scheme = c("legacy", "paper")) {
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
weight_composition <- function(x, weight_scheme = "legacy") {
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
