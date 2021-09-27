#' Load a Matrix corresponds to metacommunity composition
#'
#' This function loads a file as a matrix.
#' It assumes that columns correspond to local communities (plots or patches)
#' while the columns correspond to species.
#' The matrix entries correspond to species property (e.g, presence/absence).
#' @importFrom magrittr `%>%`
#' @import dplyr
#' @import tidyr
#' @import geometry
#' @param meta_composition metacommunity composition
#' @param weights whether the hyper points shouldbe weightedß
#' @return the geometric measure of beta_diversity
#' @export
beta_volume <- function(meta_composition, weights = T) {
  weight_composition <- function(meta_composition) {
    weight <- meta_composition %>%
      {suppressMessages(as_tibble(., .name_repair = "unique"))} %>%
      group_by(across()) %>%
      mutate(n = n()) %>%
      ungroup() %>%
      mutate(n = nrow(.) * n / sum(n)) %>%
      pull(n)
    for (i in 1:nrow(meta_composition)) {
      meta_composition[i, ] <- weight[i] * meta_composition[i, ]
    }
    meta_composition
  }

  P <- ncol(meta_composition)
  if (weights) meta_composition <- weight_composition(meta_composition)
  tryCatch(
    {
      P * (convhulln(meta_composition, output.options = TRUE)$vol)^(1/P)
    },
    error = function(e) {
      0
    }
  )
}


