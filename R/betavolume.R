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

  meta_composition <- process_meta_composition(meta_composition,
                                               weights,
                                               remove_unique)

  d <- ncol(meta_composition)

  if (missing(hypervolume_method)) {
    hypervolume_method <- ifelse(d > dim_threshold, "hyper_normal", "deterministic")
  }

  meta_composition <- rbind(meta_composition, rep(0, d))

  estiamte_volume(meta_composition, hypervolume_method, dimension = d)
}
