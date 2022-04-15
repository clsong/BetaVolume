#' Compute the patch-specific contribution to beta diversity
#'
#' @param meta_composition metacommunity composition
#' @param weights whether non-unique compositions should be weighted
#' @param remove_unique whether unique compositions should be removed
#' @return the value of the geometric measure of beta_diversity
#' @export
disentangle_patch_contribution <- function(meta_composition,
                                           weights = F,
                                           remove_unique = F) {
  meta_composition <- preprocess_meta_composition(
    meta_composition,
    weights,
    remove_unique
  )

  d <- ncol(meta_composition)

  vol_original <- betavolume(meta_composition,
    preprocess = F
  )

  1:nrow(meta_composition) %>%
    map_dbl(~ betavolume(meta_composition[-., ], preprocess = F)) %>%
    {
      vol_original - .
    } %>%
    {
      . / sum(.)
    }
}

if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))
