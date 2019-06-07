using DynamicalBilliards
import PyPlot:plot
export billiard_dual_stadium, plot

function billiard_inverse_stadium(l=1.0, w=1.0)
    l = convert(AbstractFloat, l)
    l, w = promote(l,w)

    o = typeof(l)(0.0)
    bw = FiniteWall([o,o], [l,o], [o, -w], false, "Bottom wall")
    tw = FiniteWall([l,w], [o,w], [o, w], false, "Top wall")
    leftc = Disk([o, w/2], w/2, "Left disk")
    rightc = Disk([l, w/2], w/2, "Right disk")

    return Billiard(bw, rightc, tw, leftc)
end

"""
    billiard_dual_stadium(l = 1.0, w = 1.0)
Returns a [DualBilliard](@ref) containing a stadium billiard (see 
[DynamicalBilliards.billiard_bunimovich](@ref)) and an inverse stadium.
"""
billiard_dual_stadium(l = 1.0, w = 1.0) =
    (billiard_stadium(l, w), billiard_inverse_stadium(l, w))


"""
    dual_inside(bds::DualBilliard, ω)
    dual_inside(bd::Billiard, ω)
Creates a `DualParticle` at the random
coordinates in `bd`
"""
function dual_inside(bd::Billiard{T}, ω::T) where T
    p1 = randominside(bd)
    return (p1, MagneticParticle(p1, ω))
end

dual_inside(bds::DB, ω) = dual_inside(bds[1], ω)


plot(bds::DB; kwargs...) = plot(bds[1]; kwargs...)
