#' Load a Matrix
#'
#' This function loads a file as a matrix. It assumes that the first column
#' contains the rownames and the subsequent columns are the sample identifiers.
#' Any rows with duplicated row names will be dropped with the first one being
#' kepted.
#'
#' @param infile Path to the input file
#' @return A matrix of the infile
#' @export
weight_composition <- function(meta_composition) {
  weight <- meta_composition %>%
    as_tibble() %>%
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

beta_volume <- function(meta_composition){
  tryCatch(
    {
      convhulln(comp.dat_binary[[i]], output.options=TRUE)$vol
    },
    error = function(e) {
      0
    }
  )
}
