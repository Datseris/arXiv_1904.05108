using DrWatson
quickactivate(@__DIR__, "MagneticBilliardsLyapunovs")
using DynamicalBilliards, PyPlot, LinearAlgebra


include("constants.jl")
include(srcdir()*"plot_perturbationgrowth.jl")
include(srcdir()*"unitcells.jl")

function highlight_ucells(ucells, logΔ, ts, ax = gca())
    sca(ax)
    xlim(xlim()...)
    ylim(ylim()...)
    α = 0.15
    randc = ["C$i" for i in (0,2,4,5,6,7,9)]
    if length(ucells[1]) == 3
        for (k,(s, f, c)) in enumerate(ucells)
            color = c == 2 ? "C1" : randc[(k-1)%length(randc) + 1]
            if c == 2
                axhspan(logΔ[s-1], logΔ[f],
                color = color, alpha = α, hatch = "X")
            else
                axhspan(logΔ[s-1], logΔ[f],
                color = color, alpha = α)
            end
        end
    else
        for (k, (s, f)) in enumerate(ucells)
            color = randc[(k-1)%length(randc) + 1]
            axhspan(logΔ[s-1], logΔ[f], color = color, alpha = α)
        end
    end
end

figure(figsize = (2figx, 10))
# %%

h, w, r = 1.0, 0.5, 1.0

bd1 = billiard_sinai(;setting = "periodic")
bd2 = billiard_mushroom()
data = []; limits = []; ucells = []


# PSB data
p = Particle(-0.39244660487597915, 1.7903641320527997, 1.9400545647056406)

ts, Gs, obst = perturbationgrowth(p, bd1, 200.0);
Δs = perturbationevolution(Gs)
logΔ = log.(norm.(Δs))
# colorcodeplot(ts, logΔ, obst)

push!(data, (ts, log.(norm.(Δs)), obst))
push!(limits, (85.0, 105.0))
push!(ucells, psb_unitcells(obst))

# MPSB data
p = MagneticParticle(0.6065358878346092, 0.1873249063204112,
0.9218787810732642, 2.0)

ts, Gs, obst = perturbationgrowth(p, bd1, 200.0);
Δs = perturbationevolution(Gs)
logΔ = log.(norm.(Δs))
# colorcodeplot(ts, logΔ, obst)

push!(data, (ts, log.(norm.(Δs)), obst))
push!(limits, (70.0, 90.0))
push!(ucells, psb_unitcells(obst))

# Mushroom data
p = Particle(-0.39244660487597915, 1.7903641320527997, 1.9400545647056406)

ts, Gs, obst = perturbationgrowth(p, bd2, 200.0);
Δs = perturbationevolution(Gs)
logΔ = log.(norm.(Δs))

push!(data, (ts, log.(norm.(Δs)), obst))
push!(limits, (130.0, 180.0))
push!(ucells, mushroom_unitcells(obst))


# %% Produce all colored perturbation growths:
for i in eachindex(data)
    (ts, logΔ, obst) = data[i]
    ax1 = subplot2grid((2, 3), (0,i-1))
    if i == 3
        obstsize, obstcolor = mushsize, mushcolor
    else
        obstsize, obstcolor = sinaisize, sinaicolor
    end

    colorcodeplot(ts, logΔ, obst, ax = ax1, obstsize = obstsize,
    obstcolor = obstcolor)
    u = ucells[i]
    highlight_ucells(u, logΔ, ts, ax1)
    ax2 = subplot2grid((2, 3), (1,i-1))
    if i == 1
        global axx = ax2
    end
    colorcodeplot(ts, logΔ, obst, ax = ax2, tlim = limits[i],
    obstsize = obstsize, obstcolor = obstcolor)
    highlight_ucells(u, logΔ, ts, ax2)

    ax2.set_xlabel("\$t\$")
    ax1.set_xticks([0, 200])
    if i == 1
        for ax in (ax1, ax2)
            ax.set_ylabel("\$\\log(|\\delta\\Gamma|)\$")
        end
    end
    # plot zoomin box
    x = ax2.get_xlim()
    y = ax2.get_ylim()
    rec = matplotlib.patches.Rectangle([x[1], y[1]], x[2]-x[1], y[2]-y[1],
          fill = false, color = "k", linewidth = 2.0, zorder = 99)
    ax1.add_artist(rec)

    con = matplotlib.patches.ConnectionPatch(
    xyA= (0.5(x[2]-x[1]) + x[1], y[2]), xyB= (0.5(x[2]-x[1]) + x[1], y[1]), coordsA="data", coordsB="data",
    axesA=ax2, axesB=ax1, color="black", lw = 1.0, zorder = 99
    )
    ax2.add_artist(con)
end



# Add text in panel (c)
t, G = data[1]
t1 = findfirst(x -> x > 90.5, t)
t2 = findfirst(x -> x > 98.4, t)
axx.annotate("\$|\\delta\\Gamma_{\\!\\!j}|\$", (t[t1], G[t1]),
xytext = (-18, -30), textcoords ="offset points", size = 30)

axx.annotate("\$|\\delta\\Gamma_j'|\$", (t[t1+1], G[t1+1]), xytext = (-18, 15),
textcoords ="offset points", size = 30)

axx.annotate("\$|\\delta\\Gamma_{\\!\\!j+1}|\$", (t[t2], G[t2]), xytext = (5, -8),
textcoords ="offset points", size = 30)
axx.text((t[t1] + t[t2])/2, G[t1]-4.5, "\$\\Delta t_j\$", size = 30)

axx.plot((t[t1], t[t2]), (G[t1], G[t1]), "k--")
axx.plot((t[t2], t[t2]), (G[t1], G[t2]), "k--")

tight_layout()
add_identifiers!()
subplots_adjust(wspace=0.12, hspace = 0.2, left = 0.06,right = 0.97, bottom = 0.1)

# savefig(papersdir()*"figures/pertgrowth.png")
