#' Preprocess the metacommunity composition to facilitate computation of hypervolume
#' @param meta_composition metacommunity composition
#' @param weights whether non-unique compositions should be weighted
#' @param remove_unique whether unique compositions should be removed
#' @return processed metacommunity composition
#' @export
preprocess_meta_composition <- function(meta_composition,
                                        weights = F,
                                        remove_unique = F) {
  weight_composition <- function(meta_composition) {
    meta_composition <- meta_composition[which(rowSums(meta_composition) != 0), ]
    colnames(meta_composition) <- paste0("x", 1:ncol(meta_composition))
    weight <- meta_composition %>%
      as_tibble() %>%
      group_by(across()) %>%
      mutate(n = n()) %>%
      ungroup() %>%
      mutate(n = nrow(meta_composition) * n / sum(n)) %>%
      pull(n)
    for (i in 1:nrow(meta_composition)) {
      meta_composition[i, ] <- weight[i] * meta_composition[i, ]
    }
    meta_composition
  }

  if (remove_unique) {
    meta_composition <- unique(meta_composition, MARGIN = 1)
    meta_composition <- unique(meta_composition, MARGIN = 2)
  }

  if (nrow(meta_composition) < ncol(meta_composition)) {
    meta_composition <- t(meta_composition)
  }

  d <- ncol(meta_composition)
  if (weights) meta_composition <- weight_composition(meta_composition)

  meta_composition
}
