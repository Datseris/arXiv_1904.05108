using Revise, DrWatson
quickactivate(@__DIR__, "MagneticBilliardsLyapunovs")
using InverseStadium, DynamicalBilliards, PyPlot, LinearAlgebra

bds = billiard_dual_stadium()
ω = 1.0
#pss = [dual_inside(bds, ω) for i ∈ 1:200]

# this particle should be in the chaotic sea for ω ∈ [1.0, 1000.0]
pos, vel = from_bcoords(2, -0.5, bds[1])
p1 = Particle(pos, vel)
ps = (p1, MagneticParticle(p1, ω))



Nboxes = 300
δξ = totallength(bds[1])/Nboxes
δφ = 2/Nboxes

r, dict = boundarymap_portion(bds, 2e8, ps, δξ, δφ)
psr = phasespace_portion(bds, ps, dict, δξ, δφ)

plot_boundarymap_portion(dict, δξ, δφ)

bmap = boundarymap(ps, bds, 2e8)
plot_boundarymap(bmap...; ms = 0.1, mew = 0)
