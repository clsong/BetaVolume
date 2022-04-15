#' Calculate how much the original hypervolume is changed by the similarity matrix
#'
#' @param similarity_matrix a matrix denote the level of species similarity or functional complementarity
#' @return the size of the similarity angle
#' @export
get_similarity_angle <- function(similarity_matrix){
  S <- nrow(similarity_matrix)
  omega <- function(S, Sigma) {
    m <- matrix(0, S, 1)
    a <- matrix(0, S, 1)
    b <- matrix(Inf, S, 1)
    d <- pmvnorm(lower = rep(0, S), upper = rep(Inf, S), mean = rep(0, S), sigma = Sigma)
    out <- d[1]^(1 / S)
    return(out)
  }
  f <- function(m) class(try(solve(t(m) %*% m), silent = T)) == "matrix"
  if (f(similarity_matrix) == FALSE) {
    return(0)
  }
  else {
    Sigma <- solve(t(similarity_matrix) %*% similarity_matrix)
    return(omega(S, Sigma))
  }
}
