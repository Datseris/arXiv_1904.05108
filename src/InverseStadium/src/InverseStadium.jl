module InverseStadium
using DynamicalBilliards


DualParticle{T} = NTuple{2, AbstractParticle{T}}
DualBilliard{T} = NTuple{2, Billiard{T}}

function DualParticle(p::Particle{T}, ω::T) where {T}
    return (p, MagneticParticle(p, ω))
end

const DP = DualParticle
const DB = DualBilliard

export DualParticle, DP,
    DualBilliard, DB


include("propagation.jl")
include("billiards.jl")
include("boundarymaps.jl")
include("lyapunovs.jl")

end # module
