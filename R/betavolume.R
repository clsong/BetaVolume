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
#' @import MVNH
#' @import purrr
#' @param meta_composition metacommunity composition
#' @param weights whether the hyper points shouldbe weightedß
#' @return the geometric measure of beta_diversity
#' @export
beta_volume <- function(meta_composition, weights = T, method, dim_threshold = 10) {
  if (nrow(meta_composition) < ncol(meta_composition)) {
    meta_composition <- t(meta_composition)
  }

  P <- ncol(meta_composition)
  if (weights) meta_composition <- weight_composition(meta_composition)

  meta_composition <- rbind(meta_composition, rep(0, P))

  estiamte_volume <- function(meta_composition, method, dimension){
    if(method == 'deterministic'){
      hypervolume <-tryCatch(
        {
           P *(convhulln(meta_composition, output.options = TRUE)$vol)^(1 / P)
        },
        error = function(e) {
           0
        }
      )
    }
    if(method == 'hyper_normal'){
      hypervolume <- P*(MVNH_det(meta_composition)[1])^(1/P)*4
    }
    hypervolume
  }

  if(missing(method)){
     method <- ifelse(P > dim_threshold, 'hyper_normal', 'deterministic')
  }

  estiamte_volume(meta_composition, method, dimension = P)
}

#' @export
weight_composition <- function(meta_composition) {
  meta_composition <- meta_composition[which(rowSums(meta_composition) != 0), ]
  weight <- meta_composition %>%
    {
      suppressMessages(as_tibble(., .name_repair = "unique"))
    } %>%
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

#'
#' calcualte_max_MVHH_det <- function(P){
#'   # matrix(rep(c(1,0), P), ncol = P) %>%
#'   #   split(rep(1:ncol(.), each = nrow(.))) %>%
#'   #   cross() %>%
#'   #   map(unlist) %>%
#'   #   bind_rows() %>%
#'   #   rbind(rep(0, P)) %>%
#'   #   MVNH_det() %>%
#'   #   first()
#'
#'   expand_grid(
#'     a = 0:1,
#'     b = 0:1
#'   ) %>%
#'     {.[rep(seq_len(nrow(.)), each = 2^P/4), ]} %>%
#'     cov()
#' }

# max_beta <- function(P) {
#   matrix(rep(c(1, 0), P), ncol = P) %>%
#     split(rep(1:ncol(.),
#               each = nrow(.)
#     )) %>%
#     cross() %>%
#     map(unlist) %>%
#     bind_rows()
# }
#
# max_beta(10) %>%
#   select(`1`, `2`) %>%
#   distinct() %>%
#   {.[rep(seq_len(nrow(.)), each = 256), ]} %>%
#   cov()
#
# cov(max_beta(10)) %>% diag()
#
# MVNH_det(max_beta(10))[1]
#
# MVNH_det(cov = cov(max_beta(10)), cov.matrix=TRUE)

# Maximum beta diversity
# max_MVHH_det <- 1:20 %>%
#   pro_map_dbl(calcualte_max_MVHH_det)
# max_MVHH_det
# max_MVHH_det_2 <- 21:30 %>%
#     future_map_dbl(calcualte_max_MVHH_det, .progress = T)
# max_MVHH_det <- 1:100 %>%
#   purrr::map_dbl(~MVNH_det(cov = diag(1/4, .), cov.matrix=TRUE)[1])
# usethis::use_data(max_MVHH_det, overwrite = TRUE)
