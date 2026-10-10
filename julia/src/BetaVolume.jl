"""
    BetaVolume

Geometric (hypervolume) measures of beta diversity (Song et al. 2025, *Ecological
Monographs* 95:e70008, doi:10.1002/ecm.70008). Julia port of the R package
`BetaVolume` (version 0.2.0).

A metacommunity matrix (communities × species) is embedded as a cloud of points in a
`d`-dimensional space, `d = min(#species, #communities)`. Together with the origin the
points span a convex polytope `P` and

    β_vol = d ⋅ vol(P)^(1/d).

Main functions (R names in parentheses, available as aliases):

- [`betavolume`](@ref) (`betavolume`): geometric beta diversity, exact (convex hull) or
  hyper-normal (ellipsoid) estimator.
- [`preprocess`](@ref) (`preprocess_meta_composition`): orientation and weighted
  embedding of duplicated presence/absence compositions.
- [`contributions`](@ref) (`disentangle_individual_contribution`): contribution of each
  community (or species) to beta diversity.
- [`compressed_angle`](@ref) (`get_compressed_angle`): rescaling by species similarity or
  functional complementarity.
- [`temporal_overlap`](@ref) (`beta_temporal_overlap`): temporal change of beta diversity
  as (oriented) hypervolume overlap.
"""
module BetaVolume

using LinearAlgebra
using SparseArrays
using Statistics
using Random
import DirectQhull
import MvNormalCDF

export preprocess, betavolume, contributions, compressed_angle, temporal_overlap
export preprocess_meta_composition, disentangle_individual_contribution,
       get_compressed_angle, beta_temporal_overlap

include("utils.jl")
include("qhull.jl")
include("lp.jl")
include("preprocess.jl")
include("volume.jl")
include("contributions.jl")
include("angle.jl")
include("temporal.jl")

function __init__()
    _init_qhull_errfile()
    return nothing
end

end # module
