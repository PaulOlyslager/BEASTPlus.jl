module BEASTPlus
using BEAST
using SparseArrays
using CompScienceMeshes
# Write your package code here.
include("quadstrat.jl")
include("identitymatrix.jl")
include("kernels.jl")
include("visualize_basis.jl")
include("tracespace.jl")
include("potentials/potentials.jl")
include("potentials/math_expressions.jl")
include("potentials/regularization.jl")
end
