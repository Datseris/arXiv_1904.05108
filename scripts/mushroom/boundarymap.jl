using DrWatson
quickactivate(@__DIR__, "MagneticBilliardsLyapunovs")
using DynamicalBilliards, PyPlot

# %%
h, w, r = 0.5, 0.5, 1.0
bd = billiard_mushroom(h, w, r)

n = 2000 # how many particles to create
t = 500 # how long to evolve each one

ps = [randominside(bd) for i = 1:n]
colors = [MushroomTools.is_regular(p,h,w,r) ? "#233B43" : "#E84646" for p in ps]

bmap, arcs = parallelize(boundarymap, bd, t, ps)

plot_boundarymap(bmap, arcs, color = colors, ms = 0.25)
tight_layout()
