#### extra mathematical support for potentials
import Base: +,*,-
import BEAST: PotentialOperator, potential
struct LinearCombinationOfPotentials <: PotentialOperator
    potentials::Vector{PotentialOperator}
    coefficients::Vector{}
end

+(p1::PotentialOperator,p2::PotentialOperator) = LinearCombinationOfPotentials([p1,p2],[1,1])
*(c::Number,p::PotentialOperator) = LinearCombinationOfPotentials([p],[c])
*(p::PotentialOperator,c::Number) = LinearCombinationOfPotentials([p],[c])
-(p1::PotentialOperator,p2::PotentialOperator) = LinearCombinationOfPotentials([p1,p2],[1,-1])
+(p1::LinearCombinationOfPotentials,p2::PotentialOperator) = LinearCombinationOfPotentials([p1.potentials...,p2],[p1.coefficients...,1])
+(p1::PotentialOperator,p2::LinearCombinationOfPotentials) = LinearCombinationOfPotentials([p1,p2.potentials...],[1,p2.coefficients...])
+(p1::LinearCombinationOfPotentials,p2::LinearCombinationOfPotentials) = LinearCombinationOfPotentials([p1.potentials...,p2.potentials...],[p1.coefficients...,p2.coefficients...])
*(c::Number,p::LinearCombinationOfPotentials) = LinearCombinationOfPotentials(p.potentials,c*p.coefficients)
*(p::LinearCombinationOfPotentials,c::Number) = LinearCombinationOfPotentials(p.potentials,c*p.coefficients)
-(p1::LinearCombinationOfPotentials,p2::PotentialOperator) = LinearCombinationOfPotentials([p1.potentials...,p2],[p1.coefficients...,-1])
-(p1::PotentialOperator,p2::LinearCombinationOfPotentials) = LinearCombinationOfPotentials([p1,p2.potentials...],[1,-p2.coefficients...])
-(p1::LinearCombinationOfPotentials,p2::LinearCombinationOfPotentials) = LinearCombinationOfPotentials([p1.potentials...,p2.potentials...],[p1.coefficients...,-p2.coefficients...])

function potential(op::LinearCombinationOfPotentials,points, coeffs, basis::BEAST.Space; 
    type=SVector{3,ComplexF64},
	quadstrat=nothing)
    if quadstrat === nothing
        return sum(ci*potential(opi,points,coeffs,basis,type=type) for (opi,ci) in zip(op.potentials,op.coefficients))
    else
        return sum(ci*potential(opi,points,coeffs,basis,type=type,quadstrat=quadstrat) for (opi,ci) in zip(op.potentials,op.coefficients))
    end
end
function potential(op::LinearCombinationOfPotentials,points, coeffs, basis::BEAST.DirectProductSpace; 
    type=SVector{3,ComplexF64},
	quadstrat=nothing)
    if quadstrat === nothing
        return sum(ci*potential(opi,points,coeffs,basis,type=type) for (opi,ci) in zip(op.potentials,op.coefficients))
    else
        return sum(ci*potential(opi,points,coeffs,basis,type=type,quadstrat=quadstrat) for (opi,ci) in zip(op.potentials,op.coefficients))
    end
end
