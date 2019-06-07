using DrWatson
quickactivate(@__DIR__, "MagneticBilliardsLyapunovs")
using DynamicalBilliards, PyPlot, LinearAlgebra
include("constants.jl")

figure(figsize = (figx, 6))
ax1 = subplot(1,2,1)
ax2 = subplot(1,2,2)
# %%
ax1.clear()
l, w, r = 0.5, 1.0, 1.0

bd = billiard_mushroom(l, w, r, door = false)
plot(bd, ax = ax1)
# Add background color to regions:
wedge = matplotlib.patches.Wedge([0.0, l], r;
theta1 = 0, theta2 = 180, color = "C1", alpha = 0.1)
ax1.add_artist(wedge)

rect = matplotlib.patches.Rectangle((-w/2, 0), w, l, color = "C2", alpha = 0.1)
ax1.add_artist(rect)

# plot center of billiard
ax1.scatter([0.0], [l], color = "black", s = 50)
# ax1.axis("off")
# Plot chaotic and regualr:
# p = MushroomTools.randomregular(l, w, r)
p = Particle( 0.5348673962907543, 1.3449360143733662, -1.5095927214845015)
δ = SVector(0.0, l) .- p.pos
vperp = [-p.vel[2],p.vel[1]]
ra = abs(dot(vperp, δ))

@assert MushroomTools.is_regular(p, l, w, r)
bounce!(p, bd) # be sure to start from obstacle
x, y = timeseries!(p, bd, 4)

ax1.add_artist(PyPlot.matplotlib.patches.Arc(
[0.0, l], 2ra, 2ra, theta1 = 0.0, theta2 = 180.0,
edgecolor = "C0", lw = 2.0,
fill = false, linestyle = "dashed"))

ax1.plot(x, y, color = "C0")

# Chaotic
p = Particle(0.04952884277931302,
0.8665808844454129, -0.5909640355730128)
δ = SVector(0.0, l) .- p.pos
vperp = [-p.vel[2],p.vel[1]]
ra = abs(dot(vperp, δ))

ax1.add_artist(PyPlot.plt.Circle([0.0, l], ra, edgecolor = "C3", lw = 2.0,
fill = false, linestyle = "dashed"))

@assert !MushroomTools.is_regular(p, l, w, r)

bounce!(p, bd) # be sure to start from obstacle
x, y = timeseries!(p, bd, 4)
ax1.plot(x, y, color = "C3")
ax1.set_ylim(-0.1, l+r+0.05)
ax1.set_xlim(-w-0.05, w+0.05)
ax1.axis("off")

nice_arrow(0, -0.1, w, 0, ax1; tex = "\$w\$", xo = w/2, yo = -0.05)
nice_arrow(-w/2-0.1, l/2, 0, l, ax1; tex = "\$h\$", xo = -0.15, yo = -0.05)

ax1.text(w/2 + 0.1, l/2, "stem")
ax1.text(w/2 + 0.1, l + r - 0.15, "cap")

# Add radius = 1
phi = 80 * π/180
ax1.arrow(0, l, 0.9cos(phi), 0.9sin(phi); width = 0.02,
edgecolor = "k", facecolor = "k")
ax1.text(1.01cos(phi), l+ 1.01sin(phi), "\$r=1\$")

# %% MPSB
ax2.clear()
ax2.set_aspect("equal")
bd = billiard_sinai(;setting = "periodic")

p = Particle(0.8, 0.1, -π/4)
bounce!(p, bd)
x, y = timeseries!(p, bd, 11)
plot(bd, -1, -3, 6, 3, ax = ax2)
ax2.plot(y, x, color = "C1")
# plot(p, ax = ax2)

mp = MagneticParticle(0.10368768114405014, 0.980396523525801,
-1.8830840289200883, 2.0)

@assert DynamicalBilliards.distance(mp, bd) > 0
@assert !ispinned(mp, bd)
bounce!(mp, bd)
x, y = timeseries!(mp, bd, 8)
ax2.plot(y, x, color = "C3")
ax2.set_ylim(-3, 2.2)

ax2.add_artist(PyPlot.plt.Circle([2.5, 1.5], 0.5, edgecolor = "C0", lw = 2.0,
fill = false))
ax2.add_artist(PyPlot.plt.Circle([4, 1.0], 1/0.88, edgecolor = "C4", lw = 2.0,
fill = false))
ax2.axis("off")
ax2.text(0.36, 1.4, "\$r\$", size = 32)


# %% tight layouting and axis numbers
tight_layout()
add_identifiers!()
subplots_adjust(wspace=0.15, hspace = 0.00, left = 0.05,right = 0.97, bottom = 0.00, top = 1.0)

savefig(papersdir()*"figures/billiards.png", transparent = true)
