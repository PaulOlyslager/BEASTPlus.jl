abstract type _TupleFunction end

struct Times <: _TupleFunction end
struct Cross <: _TupleFunction end
struct RealDot <: _TupleFunction end

function BEAST.potential(op::BEAST.PotentialIntegralOperator{D,U, <: _TupleFunction,W}, points, coeffs, basis; 
    type=SVector{3,ComplexF64},
	quadstrat=defaultquadstrat(op, basis)) where {D,U,W}

    return potential(PotentialIntegralOperatorKern(op.kernel,operator(op.op2)),points,coeffs,op.bfunc(basis);type,quadstrat)
end
operator(op::Times) = *
operator(op::Cross) = ×
operator(op::RealDot) = (x,y) -> transpose(x)*y