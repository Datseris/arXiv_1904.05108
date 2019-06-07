using DrWatson
include("constants.jl")
quickactivate(@__DIR__, "MagneticBilliardsLyapunovs")
using DynamicalBilliards, Statistics, BSON, LsqFit
using Measurements: measurement, value, uncertainty
import DynamicalBilliards: SV

file = datadir()*"/mushrooms/unit_cell_times.bson"
BSON.@load file result hs ws

# figure 1 with p.s. portions and trapping
figure(figsize = (figx, 12))
axtrap = subplot(2,2,1)
axkappa = subplot(2,2,2)
axgc1 = subplot(2,2,3)
axgc2 = subplot(2,2,4)

# %%
l = 1.0; w = 0.5; r = 1.0
bd = billiard_mushroom(l, w, r; door = false)
axtrap.clear()
plot(bd, ax = axtrap)
axtrap.set_ylim(-0.1, l+r+0.05)
axtrap.set_xlim(-r-0.05, r+0.05)
axtrap.axis("off")

# trapped particle
p = Particle(
    SV(0.0581140949904988, 0.0),
    SV(-0.26306015731410987, +0.9647794326341518),
    SV(0.0, 0.0)
)
x,y,vx,vy = timeseries(p, bd, 33)
ran = 2:11
raafter = 11:14
axtrap.scatter(x[ran[1]], y[ran[1]], color = "C1", zorder = 99, s = 100.0)
axtrap.scatter(x[ran[end]], y[ran[end]], color = "C1", zorder = 99, s = 100.0,
facecolors = "None", edgecolors="C1", linewidths = 4.0)
axtrap.plot(x[ran], y[ran], color = "C1", ls = "dashed", lw = 2.0)

axtrap.scatter(x[raafter[1]], y[raafter[1]], color = "C9", zorder = 99, s = 50.0)
axtrap.scatter(x[raafter[end]], y[raafter[end]], color = "C9", zorder = 99, s = 100.0,
facecolors = "None", edgecolors="C9", linewidths = 2.0)
axtrap.plot(x[raafter], y[raafter], color = "C9", lw = 3.0)
# standard chaotic
p = Particle(-0.22, 0.0, 1.2)
x,y,vx,vy = timeseries(p, bd, 20)
ran = 2:8
axtrap.plot(x[ran], y[ran], color = "C4", lw = 2.0, alpha = 0.8)
axtrap.scatter(x[ran[1]], y[ran[1]], color = "C4", zorder = 99, s = 100.0)
axtrap.scatter(x[ran[end]], y[ran[end]], color = "C4", zorder = 99, s = 100.0,
facecolors = "None", edgecolors="C4", linewidths = 2.0)
#another standard chaotic
# p = Particle(0.2, 0.0, π - π/4)
# x,y,vx,vy = timeseries(p, bd, 100)
# ran = 63:65
# axtrap.plot(x[ran], y[ran], color = "C9", lw = 2.0, alpha = 0.8)
# axtrap.scatter(x[ran[1]], y[ran[1]], color = "C9", zorder = 99, s = 100.0)
# axtrap.scatter(x[ran[end]], y[ran[end]], color = "C9", zorder = 99, s = 100.0,
# facecolors = "None", edgecolors="C9", linewidths = 2.0)

# highlight cap head:
d = bd[4]
theta1 = atan(d.facedir[2], d.facedir[1])*180/π + 90
theta2 = theta1 + 180
s1 = matplotlib.patches.Arc(d.c, 2d.r, 2d.r, theta1 = theta1, theta2 = theta2,
edgecolor = (0, 0.2, 0), lw = 2.0)
axtrap.add_artist(s1)

# add width and height
nice_arrow(0, -0.1, w, 0, axtrap; tex = "\$w\$", xo = w/2, yo = -0.05)
nice_arrow(-w/2-0.1, l/2, 0, l, axtrap; tex = "\$h\$", xo = -0.15, yo = -0.05)


# %% Stem return time
axkappa.clear()
hi = 10; hhh = hs[hi];
tc, tl = result[hi, :, 1], result[hi, :, 2]
tlc, fl = result[hi, :, 3], result[hi, :, 4]
καππα(h, w) =  (π/2*(2h + π/2)) * 2MushroomTools.V_3D_cha(h, w, 1.0) / (2π*w*(2h + π/2))
καππα2(h, w) = MushroomTools.V_3D_cha(h, w, 1.0) * 1/(2w)


axkappa.plot(ws, καππα2.(hhh, ws),
label = "analytic", color = coolcolors[1], zorder = 1,ls="dashed")
axkappa.plot(ws, @.((1-fl)*tc + fl*(tl + tlc)),
label = "numeric", color = coolcolors[2],  zorder = 99)
axkappa.legend()
axkappa.set_xticks([0, 2])
axkappa.set_xlabel("\$w\$", labelpad = -10)
axkappa.set_ylabel("\$\\kappa\$", labelpad = -10)
axkappa.set_yticks(5.6:0.3:6.5)

nice_arrow(0, -0.1, w, 0, axtrap; tex = "\$w\$", xo = w/2, yo = -0.05)
nice_arrow(-w/2-0.1, l/2, 0, l, axtrap; tex = "\$h\$", xo = -0.15, yo = -0.05)

# %% phase space portions
axgc1.clear()
axgc2.clear()
let
    ws = 0.0:0.01:2.0
    for (j, h) in enumerate([0.5, 1.0, 1.5])
        gc = MushroomTools.V_3D_cha.(h, ws, 1.0)
        axgc2.plot(ws, gc, label = "\$h = $h\$", color = coolcolors[j])
    end
    hs = 0.0:0.01:2.0
    for (j, w) in enumerate([0.5, 1.0, 1.5])
        gc = MushroomTools.V_3D_cha.(hs, w, 1.0)
        axgc1.plot(ws, gc, label = "\$w = $w\$", color = coolcolors[j])
    end
    for ax in (axgc1, axgc2)
        ax.set_ylabel("\$V_\\mathrm{CH}\$")
        ax.legend()
    end
end

axgc1.set_xlabel("\$h\$", labelpad = -10)
axgc2.set_xlabel("\$w\$", labelpad = -10)
for ax in (axgc1, axgc2)
    ax.set_yticks(0:10:30)
    ax.set_xticks([0,2])
    ax.set_ylim(0,30)
end

# %%
tight_layout()
add_identifiers!()
subplots_adjust(wspace=0.3, hspace = 0.25, left = 0.1,right = 0.97, bottom = 0.08)
savefig(papersdir()*"figures/mushroom1.png", transparent = true)
