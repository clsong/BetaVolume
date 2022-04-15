#' Load a Matrix corresponds to metacommunity composition
#'
#' This function loads a file as a matrix.
#' It assumes that columns correspond to local communities (plots or patches)
#' while the columns correspond to species.
#' The matrix entries correspond to species property (e.g, presence/absence).
#'
#' @param meta_composition metacommunity composition
#' @param weights whether non-unique compositions should be weighted
#' @param hypervolume_method which method to estimate the hypervolume
#' @param dim_threshold the threshold of which default hypervolume_method is used
#' @param remove_unique whether unique compositions should be removed
#' @return the value of the geometric measure of beta_diversity
#' @export
betavolume <- function(meta_composition,
                       weights = F,
                       hypervolume_method,
                       dim_threshold = 10,
                       remove_unique = F) {
  weight_composition <- function(meta_composition) {
    meta_composition <- meta_composition[which(rowSums(meta_composition) != 0), ]
    weight <- meta_composition %>%
      as_tibble(.name_repair = "unique") %>%
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

  estiamte_volume <- function(meta_composition, hypervolume_method, dimension) {
    if (hypervolume_method == "deterministic") {
      hypervolume <- tryCatch(
        {
          d * (convhulln(meta_composition, output.options = TRUE)$vol)^(1 / d)
        },
        error = function(e) {
          0
        }
      )
    }
    if (hypervolume_method == "hyper_normal") {
      hypervolume_raw <- eigen(cov(meta_composition), only.values = T)$values %>%
        log() %>%
        {
          sum(.data) / length(.data)
        } %>%
        exp()

      hypervolume <- d * hypervolume_raw * 4
    }
    hypervolume
  }

  if(remove_unique){
    meta_composition <- unique(meta_composition, MARGIN = 1)
    meta_composition <- unique(meta_composition, MARGIN = 2)
  }

  if (nrow(meta_composition) < ncol(meta_composition)) {
    meta_composition <- t(meta_composition)
  }

  d <- ncol(meta_composition)
  if (weights) meta_composition <- weight_composition(meta_composition)

  if (missing(hypervolume_method)) {
    hypervolume_method <- ifelse(d > dim_threshold, "hyper_normal", "deterministic")
  }

  meta_composition <- rbind(meta_composition, rep(0, d))

  estiamte_volume(meta_composition, hypervolume_method, dimension = d)
}
