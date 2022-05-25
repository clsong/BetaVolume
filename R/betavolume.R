#' Load a Matrix corresponds to metacommunity composition
#'
#' This function loads a file as a matrix.
#' It assumes that columns correspond to local communities (plots or patches)
#' while the columns correspond to species.
#' The matrix entries correspond to species property (e.g, presence/absence).
#'
#' @param meta_composition metacommunity composition
#' @param preprocess Whether weights and remove_unique should be executed
#' @param weights whether non-unique compositions should be weighted
#' @param hypervolume_method which method to estimate the hypervolume
#' @param dim_threshold the threshold of which default hypervolume_method is used
#' @param remove_unique whether unique compositions should be removed
#' @param compressed_angle how much the original hypervolume is changed by the similarity matrix
#' @param sparsity whether the meta_composition is encoded as a sparse matrix
#' @return the value of the geometric measure of beta_diversity
#' @export
betavolume <- function(meta_composition,
                       preprocess = T,
                       weights = F,
                       remove_unique = F,
                       hypervolume_method,
                       compressed_angle = 1,
                       dim_threshold = 10,
                       sparsity = F) {
  estiamte_volume <- function(meta_composition, hypervolume_method, dimension, sparsity) {
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
    if (hypervolume_method == "hyper_normal" & !sparsity) {
      hypervolume_raw <- eigen(cov(meta_composition), only.values = T)$values %>%
        log() %>%
        {
          sum(.) / length(.)
        } %>%
        exp()

      hypervolume <- d * hypervolume_raw * 4
    }
    if (hypervolume_method == "hyper_normal" & sparsity) {
      sparse.cov <- function(x) {
        n <- nrow(x)
        m <- ncol(x)
        ii <- unique(x@i) + 1 # rows with a non-zero element

        Ex <- colMeans(x)
        nozero <- as.vector(x[ii, ]) - rep(Ex, each = length(ii)) # colmeans

        covmat <- (crossprod(matrix(nozero, ncol = m)) +
          crossprod(t(Ex)) * (n - length(ii))
        ) / (n - 1)
        covmat
      }
      cov_mat <- sparse.cov(meta_composition)
      hypervolume_raw <- eigen(cov_mat, only.values = T)$values %>%
        log() %>%
        {
          sum(.) / length(.)
        } %>%
        exp()
      d * hypervolume_raw * 4
    }
    hypervolume
  }

  if (preprocess) {
    meta_composition <- preprocess_meta_composition(
      meta_composition,
      weights,
      remove_unique
    )
  }

  d <- ncol(meta_composition)

  if (missing(hypervolume_method)) {
    hypervolume_method <- ifelse(d > dim_threshold, "hyper_normal", "deterministic")
  }

  meta_composition <- rbind(meta_composition, rep(0, d))

  compressed_angle * estiamte_volume(meta_composition, hypervolume_method, dimension = d, sparsity)
}
