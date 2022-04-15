#' Load a Matrix corresponds to metacommunity composition
#'
#' This function loads a file as a matrix.
#' It assumes that columns correspond to local communities (plots or patches)
#' while the columns correspond to species.
#' The matrix entries correspond to species property (e.g, presence/absence).
#'
#' @param meta_composition metacommunity composition
#' @param weights whether non-unique compositions should be weighted
#' @param method which method to use to compute the hypervolume
#' @param dim_threshold the threshold of which default method is used
#'
#' @importFrom geometry convhulln
#' @importFrom stats cov
#' @importFrom tibble as_tibble
#'
#' @return the value of the geometric measure of beta_diversity
#' @export
beta_volume <- function(meta_composition,
                        weights = T,
                        method,
                        dim_threshold = 10) {
  if (nrow(meta_composition) < ncol(meta_composition)) {
    meta_composition <- t(meta_composition)
  }

  P <- ncol(meta_composition)
  if (weights) meta_composition <- weight_composition(meta_composition)

  meta_composition <- rbind(meta_composition, rep(0, P))

  estiamte_volume <- function(meta_composition, method, dimension) {
    if (method == "deterministic") {
      hypervolume <- tryCatch(
        {
          P * (convhulln(meta_composition, output.options = TRUE)$vol)^(1 / P)
        },
        error = function(e) {
          0
        }
      )
    }
    if (method == "hyper_normal") {
      meta_composition <- unique(meta_composition, MARGIN = 1)
      meta_composition <- unique(meta_composition, MARGIN = 2)

      if (nrow(meta_composition) < ncol(meta_composition)) {
        meta_composition <- t(meta_composition)
      }

      P <- ncol(meta_composition)

      hypervolume_raw <- eigen(cov(meta_composition), only.values= T)$values %>%
        log() %>%
        {sum(.)/length(.)} %>%
        exp()

      hypervolume <- P * hypervolume_raw * 4
    }
    hypervolume
  }

  if (missing(method)) {
    method <- ifelse(P > dim_threshold, "hyper_normal", "deterministic")
  }

  estiamte_volume(meta_composition, method, dimension = P)
}

#' @export
weight_composition <- function(meta_composition) {
  meta_composition <- meta_composition[which(rowSums(meta_composition) != 0), ]
  weight <- meta_composition %>%
    {
      suppressMessages(as_tibble(.data, .name_repair = "unique"))
    } %>%
    group_by(across()) %>%
    mutate(n = n()) %>%
    ungroup() %>%
    mutate(n = nrow(.data) * n / sum(n)) %>%
    pull(n)
  for (i in 1:nrow(meta_composition)) {
    meta_composition[i, ] <- weight[i] * meta_composition[i, ]
  }
  meta_composition
}
