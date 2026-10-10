#' BetaVolume: geometric (hypervolume) measures of beta diversity
#'
#' @description
#' BetaVolume implements the geometric approach to beta diversity of Song et
#' al. (2025). A metacommunity matrix (communities x species) is embedded as a
#' cloud of points in a hyper-dimensional space, and beta diversity is the
#' rescaled hypervolume of the geometric object spanned by the points and the
#' origin. In contrast with classic measures, adding a distinct, species-rich
#' community to a metacommunity increases beta diversity, so the measure
#' explicitly accounts for species associations and non-additivity.
#'
#' @section Main functions:
#' \describe{
#'   \item{[betavolume()]}{Geometric beta diversity \eqn{\beta_{vol}}, exact
#'     (convex hull) or approximated (hyper-normal ellipsoid).}
#'   \item{[preprocess_meta_composition()]}{Orientation of the matrix and
#'     weighted embedding of duplicated presence/absence compositions.}
#'   \item{[disentangle_individual_contribution()]}{Contribution of each
#'     community (or species) to beta diversity.}
#'   \item{[get_compressed_angle()]}{Rescaling of the hypervolume by species
#'     similarity or functional complementarity.}
#'   \item{[beta_temporal_overlap()]}{Temporal change of beta diversity as
#'     (oriented) hypervolume overlap.}
#' }
#'
#' See `vignette("BetaVolume")` for a guided tour that reproduces the worked
#' examples of the paper.
#'
#' @references Song, C., Lu, M., Bennett, J. R., Gilbert, B., Fortin, M.-J., &
#'   Gonzalez, A. (2025). A geometric approach to beta diversity. *Ecological
#'   Monographs*, 95(1), e70008. \doi{10.1002/ecm.70008}
#'
#' @keywords internal
"_PACKAGE"
