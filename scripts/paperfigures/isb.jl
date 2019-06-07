using DrWatson
quickactivate(@__DIR__, "MagneticBilliardsLyapunovs")
include("constants.jl")
using DynamicalBilliards, InverseStadium, PyPlot
using DynamicalBilliards: cossin

#
figure(figsize = (figx, 12))
axorbit = subplot(2,1,1)
axbmap = subplot(2,2,3)
axλ = subplot(2,2,4)
l = 0.5; w = 0.5
ω = 10.0
bds = billiard_dual_stadium(l, w)
pss = [dual_inside(bds, ω) for i ∈ 1:300]
bmap = boundarymap(pss, bds, 5000.0)

# %% Plot billiard
axorbit.clear()
plot(bds[1]; ax = axorbit)
ps = dual_inside(bds, ω)
ps[1].pos = ps[2].pos = (l/2, w/2)
ps[1].vel = ps[2].vel = cossin.(π/6)
plot(ps[1])
xt, yt = timeseries!(ps, bds, 3.0)
plot(ps[1])
axorbit.plot(xt, yt)
axorbit.set_xlim(-w, l+w)
axorbit.set_ylim(-w/2, 3w/2)
axorbit.axis("off")
cir = matplotlib.patches.Circle((-0.05, 0.15), 1/ω, fill = false, color = "k", lw = 2.0)
axorbit.add_artist(cir)
axorbit.text(-0.05, 0.18, "\$1/\\omega\$", size = 16, ha = "center", va = "center")
axorbit.arrow(-0.05, 0.15, 1/ω - 0.02, 0, width = 0.005, color = "k")

# Add indication for ξ :
arc = matplotlib.patches.Arc((0, l/2), 1/(4w) + 0.05, 1/(4w) + 0.05,
theta1 = 200, theta2 = 270, ls = "dashed", lw = 2.0)
axorbit.add_artist(arc)
axorbit.arrow(0, -0.025, l/4, 0, width = 0.005, color = "k")
axorbit.text(0, -0.1, "\$\\xi\$")
# %% Plot bmap
axbmap.clear()
plot_boundarymap(bmap..., bordercolor = (0,0,0,0),
obstacleindices = false, ax = axbmap, color = "k", ms = 0.1, mew = 0.0)

# %%
# Plot lyapunovs:
data = load(datadir()*"inverse_stadium/output_Tbox=100000.0_Nbox=300.bson")
ωs = data[:ωs]
λs = data[:λss]
Vc = data[:psv_vis]

axλ.clear()
axλ.semilogx(ωs, λs ./ maximum(λs), color = coolcolors[1], label = "\$\\lambda\$")
axλ.semilogx(ωs, 1 ./ Vc ./ maximum( 1 ./ Vc), color = coolcolors[2], ls = "dashed", label = "\$1/V_C\$")
axλ.legend()
axλ.set_xlabel("\$\\omega\$")

# %% labels etc.
tight_layout()
add_identifiers!()
subplots_adjust(wspace=0.2, hspace = 0.0, left = 0.12,right = 0.97, bottom = 0.08)
savefig(papersdir()*"figures/isb.png", transparent = true)
