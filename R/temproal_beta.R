#' Compute the temporal change in beta diversity
#' @param meta_composition_1 metacommunity at time 1
#' @param meta_composition_2 metacommunity at time 2
#' @param preprocess whether to preprocess the data
#' @param ... Arguments to be passed to methods
#' @return the value of the geometric measure of beta_diversity
#' @export
beta_temporal_overlap <- function(meta_composition_1,
                                  meta_composition_2,
                                  preprocess = T,
                                  ...) {
  if (preprocess) {
    meta_composition_1 <- preprocess_meta_composition(
      meta_composition_1, ...
    )
    meta_composition_2 <- preprocess_meta_composition(
      meta_composition_2, ...
    )
  }

  d <- ncol(meta_composition_1)

  (intersectn(meta_composition_1, meta_composition_2)$ch$vol)^(1/d) / (convhulln(meta_composition_1, output.options = TRUE)$vol)^(1 / d)
}
