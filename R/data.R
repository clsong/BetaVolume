#' Example metacommunity from A. H. Gentry's forest transect data
#'
#' Presence/absence of 37 woody plant species in the 10 subplots (2 x 50 m
#' transects) of one 0.1-ha forest plot of A. H. Gentry's global forest
#' transect data set. Each plot can be viewed as a metacommunity of 10 local
#' communities; Song et al. (2025) analysed 198 such plots along a latitudinal
#' gradient.
#'
#' @format A numeric 0/1 matrix with 37 rows (species, named
#'   `FAMILY_GENUS_SPECIES`) and 10 columns (subplots `"1"` to `"10"`). A 1
#'   means that the species was recorded in the subplot.
#'
#' Because there are fewer subplots than species, the subplots define the
#' 10 axes of the embedding and the 37 species are the embedded points (see
#' [preprocess_meta_composition()]).
#'
#' @source Gentry, A. H. (1988). Changes in plant community diversity and
#'   floristic composition on environmental and geographical gradients.
#'   *Annals of the Missouri Botanical Garden*, 75, 1-34.
#'   \doi{10.2307/2399464}. Data as distributed with Kraft et al. (2011)
#'   \doi{10.1126/science.1208584} and Xing & He (2021)
#'   \doi{10.5281/zenodo.4270202}.
#'
#' @references Song, C., Lu, M., Bennett, J. R., Gilbert, B., Fortin, M.-J., &
#'   Gonzalez, A. (2025). A geometric approach to beta diversity. *Ecological
#'   Monographs*, 95(1), e70008. \doi{10.1002/ecm.70008}
#'
#' @examples
#' dim(example_meta_composition)
#' colSums(example_meta_composition) # species richness of each subplot
#' betavolume(example_meta_composition)
"example_meta_composition"
