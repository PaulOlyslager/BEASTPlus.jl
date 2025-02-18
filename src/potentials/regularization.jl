

regularize(a) = a
function regularize(a::LinearCombinationOfPotentials)
    out = regularize(a.potentials[1])*a.coefficients[1]
    for i in 2:length(a.potentials)
        out += regularize(a.potentials[i])*a.coefficients[i]
    end
    return out
end

function regularize(pot::BEAST.PotentialIntegralOperator{3,<:BEAST.HH3DGradGreen, <: Cross,W}) where {W}
    return -regularize(BEAST.PotentialOperator{2}(BEAST.HH3DGreen(pot.kernel.gamma),Times(),strace(pot.bfunc))) +
        regularize(BEAST.PotentialOperator{3}(BEAST.HH3DGreen(pot.kernel.gamma),Times(),curl(pot.bfunc)))
end

function regularize(pot::BEAST.PotentialIntegralOperator{3,<:BEAST.HH3DGreen, <: Times,W}) where {W}
    return -regularize(BEAST.PotentialOperator{2}(HH3DInt1,Times(),ntrace(pot.bfunc))) +
        regularize(BEAST.PotentialOperator{3}(HH3DInt1,Times(),divergence(pot.bfunc))) +
        regularize(BEAST.PotentialOperator{2}(HH3DInt1,Cross(),strace(pot.bfunc))) -
        regularize(BEAST.PotentialOperator{3}(HH3DInt1,Cross(),curl(pot.bfunc)))
end
function regularize(pot::BEAST.PotentialIntegralOperator{3,<: HH3DGradDivGreen, <: Times,W}) where {W}
    @warn "assumed basis function is div conforming and normal zero on boundary"
    @warn "assumed basis functio is of order 0"
    return -regularize(BEAST.PotentialOperator{2}(BEAST.HH3DGreen(pot.kernel.gamma),Times(),ntimestrace(pot.bfunc))) 
end