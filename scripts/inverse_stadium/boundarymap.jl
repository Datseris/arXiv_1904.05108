using Revise, DrWatson
quickactivate(@__DIR__, "MagneticBilliardsLyapunovs")
using InverseStadium, DynamicalBilliards, PyPlot

bds = billiard_dual_stadium()
pss = [dual_inside(bds, 10.0) for i ∈ 1:200]


bmap = boundarymap(pss, bds, 5000.0)

plot_boundarymap(bmap..., color = "k", ms = 0.1, mew = 0.0)
