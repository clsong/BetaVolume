#' Hypervolume rescaling induced by species similarity or complementarity
#'
#' @description
#' Species similarity (or functional complementarity) is a linear
#' transformation of the axes of the embedding space: with a similarity matrix
#' \eqn{S} the axis of species \eqn{j} points along the \eqn{j}-th column of
#' \eqn{S} instead of along the \eqn{j}-th unit vector (Song et al. 2025,
#' Fig. 7). The hypervolume is then rescaled by the normalised solid angle of
#' the cone spanned by the transformed axes (Eq. 8),
#' \deqn{\Omega(S) = \frac{2^d\,|\det S|}{(2\pi)^{d/2}}
#'   \int_{x \ge 0} e^{-\frac{1}{2} x^\top S^\top S x}\, dx
#'   = 2^d\; \Pr\left(Z \ge 0\right),\quad Z \sim \mathcal{N}\left(0, (S^\top S)^{-1}\right),}{
#'   Omega(S) = 2^d * P(Z >= 0), Z ~ N(0, (S'S)^-1),}
#' and beta diversity becomes
#' \eqn{\beta_{vol} = d\,(\Omega(S)\,\mathrm{vol}(P))^{1/d} = \Omega(S)^{1/d} \times d\,\mathrm{vol}(P)^{1/d}}{beta_vol = d (Omega(S) vol(P))^(1/d) = Omega(S)^(1/d) * d vol(P)^(1/d)}
#' (Eq. 9). This function returns the multiplicative factor
#' \eqn{\Omega(S)^{1/d}}, to be passed to the `compressed_angle` argument of
#' [betavolume()].
#'
#' @details
#' \itemize{
#'   \item \eqn{S = I} (totally dissimilar species): \eqn{\Omega = 1}, beta
#'     diversity is unchanged.
#'   \item Positive off-diagonal similarities shrink the hypervolume
#'     (\eqn{\Omega < 1}); if all species are totally similar (\eqn{S} a matrix
#'     of ones) the axes collapse and \eqn{\Omega = 0}.
#'   \item Negative off-diagonal entries encode functional complementarity and
#'     expand the hypervolume (\eqn{\Omega > 1}).
#' }
#' \eqn{S} need not be symmetric. Following Eq. 8, the transformed axes are the
#' columns of \eqn{S}.
#'
#' For \eqn{d \le 3}{d <= 3} the orthant probability has a closed form and the result
#' is exact. For \eqn{d > 3} it is estimated with the randomised quasi-Monte
#' Carlo algorithm of Genz and Bretz ([mvtnorm::pmvnorm()]) with a *relative*
#' error tolerance (`releps`). The orthant probability is of order \eqn{2^{-d}},
#' so the absolute tolerance used by BetaVolume 0.1.0 (0.001) was larger than
#' the quantity itself for ten or more species. Call [set.seed()] first for
#' reproducible results when \eqn{d > 3}.
#'
#' The embedding dimension of the similarity transformation is the number of
#' species, so this correction is meant for metacommunities with at least as
#' many communities as species (species are the axes of the embedding).
#'
#' @param similarity_matrix A square numeric matrix \eqn{S} whose entry
#'   \eqn{s_{ij}} gives how similar species \eqn{i} is to species \eqn{j}
#'   (1 on the diagonal; 0 = totally dissimilar, 1 = totally similar; negative
#'   values encode functional complementarity).
#' @param releps Relative error tolerance of the numerical integration (used
#'   only when the dimension exceeds 3).
#' @param maxpts Maximum number of function evaluations of the numerical
#'   integration (used only when the dimension exceeds 3).
#'
#' @return The factor \eqn{\Omega(S)^{1/d}} (a non-negative number; 1 for the
#'   identity matrix, 0 for a singular matrix). For \eqn{d > 3} the estimated
#'   absolute error of \eqn{\Omega(S)} is attached as attribute `"error"`.
#'
#' @references Song, C., Lu, M., Bennett, J. R., Gilbert, B., Fortin, M.-J., &
#'   Gonzalez, A. (2025). A geometric approach to beta diversity. *Ecological
#'   Monographs*, 95(1), e70008. \doi{10.1002/ecm.70008}
#'
#' Song, C., Rohr, R. P., & Saavedra, S. (2018). A guideline to study the
#' feasibility domain of multi-trophic and changing ecological communities.
#' *Journal of Theoretical Biology*, 450, 30-36. \doi{10.1016/j.jtbi.2018.04.030}
#'
#' Genz, A., & Bretz, F. (2009). *Computation of Multivariate Normal and t
#' Probabilities*. Springer. \doi{10.1007/978-3-642-01689-9}
#'
#' @seealso [betavolume()]
#' @export
#' @examples
#' get_compressed_angle(diag(2)) # 1: dissimilar species
#' get_compressed_angle(matrix(1, 2, 2)) # 0: identical species
#'
#' # Fig. 7C: somewhat similar species shrink the hypervolume
#' S <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
#' get_compressed_angle(S)
#'
#' # Fig. 7D: functional complementarity expands it
#' get_compressed_angle(matrix(c(1, -0.5, -0.5, 1), 2, 2))
#'
#' # Use the factor in betavolume()
#' meta <- rbind(c(1, 0), c(0, 1), c(1, 1))
#' betavolume(meta, compressed_angle = get_compressed_angle(S))
get_compressed_angle <- function(similarity_matrix, releps = 1e-4, maxpts = 1e6) {
  S <- similarity_matrix
  if (is.data.frame(S)) S <- as.matrix(S)
  if (inherits(S, "Matrix")) S <- as.matrix(S)
  if (!is.matrix(S) || !is.numeric(S) || nrow(S) != ncol(S) || nrow(S) == 0L) {
    stop("`similarity_matrix` must be a non-empty square numeric matrix.", call. = FALSE)
  }
  check_finite(S, "similarity_matrix")
  d <- nrow(S)

  # A singular transformation collapses the hypervolume.
  if (rcond(S) < sqrt(.Machine$double.eps)) {
    return(0)
  }
  if (d == 1L) {
    return(1)
  }
  omega <- solid_angle_fraction(S, releps = releps, maxpts = maxpts)
  value <- (2^d * omega)^(1 / d)
  if (!is.null(attr(omega, "error"))) {
    attr(value, "error") <- 2^d * attr(omega, "error")
  }
  value
}

#' Fraction of R^d covered by the cone spanned by the columns of S.
#'
#' Equals Pr(Z >= 0) for Z ~ N(0, (S'S)^{-1}). Closed forms for d = 2, 3;
#' Genz-Bretz quasi-Monte Carlo otherwise.
#' @noRd
solid_angle_fraction <- function(S, releps = 1e-4, maxpts = 1e6) {
  d <- nrow(S)
  sigma <- solve(crossprod(S))
  rho <- stats::cov2cor((sigma + t(sigma)) / 2)
  if (d == 2L) {
    return(1 / 4 + asin(rho[1, 2]) / (2 * pi))
  }
  if (d == 3L) {
    return(1 / 8 + (asin(rho[1, 2]) + asin(rho[1, 3]) + asin(rho[2, 3])) / (4 * pi))
  }
  p <- mvtnorm::pmvnorm(
    lower = rep(0, d), upper = rep(Inf, d), corr = rho,
    algorithm = mvtnorm::GenzBretz(maxpts = maxpts, abseps = 0, releps = releps)
  )
  out <- as.numeric(p)
  attr(out, "error") <- attr(p, "error")
  out
}
