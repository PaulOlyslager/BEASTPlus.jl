

regularize(a) = a
function regularize(a::LinearCombinationOfPotentials)
    out = regularize(a.potentials[1])*a.coefficients[1]
    for i in 2:length(a.potentials)
        out += regularize(a.potentials[i])*a.coefficients[i]
    end
    return out
end

function regularize(pot::BEAST.PotentialIntegralOperator{3,<:BEAST.HH3DGradGreen, <: Cross,W}) where {W}
    return (-1)*regularize(BEAST.PotentialIntegralOperator{2}(BEAST.HH3DGreen(pot.kernel.gamma),Times(),b->strace(pot.bfunc(b)))) +
        regularize(BEAST.PotentialIntegralOperator{3}(BEAST.HH3DGreen(pot.kernel.gamma),Times(),b->curl(pot.bfunc(b))))
end

# function regularize(pot::BEAST.PotentialIntegralOperator{3,<:BEAST.HH3DGreen, <: Times,W}) where {W}
#     return (-1)*regularize(BEAST.PotentialIntegralOperator{2}(HH3DInt1(pot.kernel.gamma),Times(),b->ntrace(pot.bfunc(b)))) +
#         regularize(BEAST.PotentialIntegralOperator{3}(HH3DInt1(pot.kernel.gamma),Times(),b->divergence(pot.bfunc(b)))) +
#         regularize(BEAST.PotentialIntegralOperator{2}(HH3DInt1(pot.kernel.gamma),Cross(),b->strace(pot.bfunc(b)))) -
#         regularize(BEAST.PotentialIntegralOperator{3}(HH3DInt1(pot.kernel.gamma),Cross(),b->curl(pot.bfunc(b))))
# end
function regularize(pot::BEAST.PotentialIntegralOperator{3,<:BEAST.HH3DGreen, <: Times,W}) where {W}
    @warn "assumed basis function is normal zero on boundary"
    @warn "assumes curl is zero of basis function"
    return regularize(BEAST.PotentialIntegralOperator{3}(HH3DInt1(pot.kernel.gamma),Times(),b->divergence(pot.bfunc(b)))) +
        regularize(BEAST.PotentialIntegralOperator{2}(HH3DInt1(pot.kernel.gamma),Cross(),b->strace(pot.bfunc(b))))
end
function regularize(pot::BEAST.PotentialIntegralOperator{3,<: HH3DGradDivGreen, <: Times,W}) where {W}
    @warn "assumed basis function is div conforming and normal zero on boundary"
    @warn "assumed basis functio is of order 0"
    return (-1)*regularize(BEAST.PotentialIntegralOperator{2}(BEAST.HH3DGreen(pot.kernel.gamma),Times(),b->ntimestrace(pot.bfunc(b)))) 
end