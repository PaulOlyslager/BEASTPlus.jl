using BEAST
using StaticArrays
const i4pi = 1/(4*pi)
struct HH3DGradGreen2{T} <: BEAST.Kernel{T}
    gamma::T
end
function (op::HH3DGradGreen2)(x,y)
    gamma = op.gamma

    r = cartesian(x) - cartesian(y)
    R = norm(r)
    gradgreen = -((expm1(-gamma*R)*(1+gamma*R)+gamma*R)/(4*pi*R^3))*r
    return gradgreen
end
struct HH3DGreen2{T} <: BEAST.Kernel{T}
    gamma::T
end
function (op::HH3DGreen2)(x,y)
    gamma = op.gamma

    r = cartesian(x) - cartesian(y)
    R = norm(r)
    iR = 1/R
    green = expm1(-gamma*R)*(i4pi*iR)
    green
end

struct HH3DGradGreenFF{T} <: BEAST.Kernel{T}
    gamma::T
end

struct HH3DGreenFF{T} <: BEAST.Kernel{T}
    gamma::T
end
struct HH3DGradGreenFF2{T} <: BEAST.Kernel{T}
    gamma::T
end

struct HH3DGreenFF2{T} <: BEAST.Kernel{T}
    gamma::T
end

struct HH3DGradDivGreen{T} <: BEAST.Kernel{T}
    gamma::T
end
function (op::HH3DGradDivGreen)(x,y)
    @error "do not integrate HH3DGradDivGreen directly"
end
struct HH3DInt1{T} <: BEAST.Kernel{T}
    gamma::T
end
struct HH3DInt2{T} <: BEAST.Kernel{T}
    gamma::T
end
# function (op::HH3DGreenFF)(x,y)
#     gamma = op.gamma

#     r = cartesian(x) - cartesian(y)
#     R = norm(r)
#     iR = 1/R
#     green = exp(-gamma*R)*(i4pi*iR)
#     green
# end
function (op::HH3DGreenFF)(x::Union{SVector,Vector},y)
    gamma = op.gamma
    x = normalize(x)
    green = exp(gamma*(x⋅cartesian(y)))*(i4pi)
    green
end
function (op::HH3DGradGreenFF)(x::Union{SVector,Vector},y)
    gamma = op.gamma
    x = normalize(x)
    green = exp(gamma*(x⋅cartesian(y)))*(i4pi)
    gradgreen = -gamma  * green * x
    gradgreen
end

function (op::HH3DGreenFF2)(x::Union{SVector,Vector},y)
    gamma = op.gamma
    x = normalize(x)
    green = expm1(gamma*(x⋅cartesian(y)))*(i4pi)
    green
end
function (op::HH3DGradGreenFF2)(x::Union{SVector,Vector},y)
    gamma = op.gamma
    x = normalize(x)
    green = expm1(gamma*(x⋅cartesian(y)))*(i4pi)
    gradgreen = -gamma  * green * x
    gradgreen
end
export HH3DGradGreen2, HH3DGreen2, HH3DGradGreenFF, HH3DGreenFF, HH3DGradGreenFF2, HH3DGreenFF2, HH3DGradDivGreen
#### exp(-ikr)/r omitted