using DrWatson
quickactivate(@__DIR__, "MagneticBilliardsLyapunovs")
include("constants.jl")
using DynamicalBilliards, Statistics, BSON, LsqFit
using Measurements: measurement, value, uncertainty

include(srcdir()*"mpsb_cellstatistics.jl")

BSON.@load datadir()*"periodic_sinai/Λ_N=5000.bson" Λ Λσ rs Bs
BSON.@load datadir()*"periodic_sinai/κ_N=5000.bson" κ κσ

ridx = (1, 3, 5)
Bs2 = 0.01:0.01:Bs[end]
Bidxs = [findfirst(isequal(B), Bs) for B in Bs2]

Bplot = collect(0.2:0.2:1.0)
Bplot = [findfirst(isequal(B), Bs2) for B in Bplot]

# Gather stats
file = datadir()*"periodic_sinai/mpsb_stats.bson"

if !isfile(file)
    stats = []
    for k in ridx
        r = rs[k]
        Δts = Vector{Float64}[];
        as = copy(Δts); δΓs = copy(Δts)
        for B in Bs2
            println(B)
            tc, inc_c, jhs = manystatistics(r, B, 2000; T = 500)
            dd = findall(isinf, exp.(inc_c))
            deleteat!(inc_c, dd)
            deleteat!(tc, dd)
            deleteat!(jhs, findall(isinf, jhs))
            push!(Δts, tc); push!(as, jhs); push!(δΓs, inc_c)
        end
        push!(stats, [Δts, as, δΓs])
    end
    BSON.@save file stats
else
    BSON.@load file stats
end

# %% Make figure
figure(figsize = (figx, 20))
axλ = subplot(3,1,1)
axa = subplot(3,2,3)
axgc = subplot(3,2,4)
axΓ = subplot(3,2,5)
zoomin = subplot(3,2,6)

function add_grid!(ax, Bs = 0:0.2:1.0; kwargs...)
    for B in Bs
        ax.axvline(B, color = "gray", alpha = 0.5, lw = 1.0)
    end
end

# %% Toy model
z_MPSB(κ,B) = sqrt(1 + ((1 - cos(2B*κ))/2B)^2 + (sin(2B*κ)/(2B))^2)
toymodel(κ, loga, B) = (loga + log(z_MPSB(κ,B)))/κ


axλ.clear()
axgc.clear()


for ax in (axλ, axgc); add_grid!(ax); end

for i in eachindex(ridx)
    Δts, as, δΓs = stats[i]

    aμ = mean.(as)
    # aσ = std.(as)
    logam = mean.([log.(z) for z in as])
    toy = toymodel.(κ[Bidxs, ridx[i]], logam, Bs2)

    axλ.plot(Bs, Λ[:, ridx[i]], color = coolcolors[i])
    # coolfill(Bs, Λ[:, ridx[i]], Λσ[:, ridx[i]], axλ, coolcolors[i], "r=$(rs[ridx[i]])")

    axλ.plot(Bs2, toy, color = coolcolors[i], ls = "dashed")
    axgc.plot(Bs, κ[:, ridx[i]] ./ maximum(κ[:, ridx[i]]), color = coolcolors[i],
    label =  "\$r=$(rs[ridx[i]])\$")
end
axλ.set_ylim(0.2, 2.5)
for ax in (axλ, axgc)
    ax.set_xticks([0, 1])
    ax.set_xlabel("\$B\$", labelpad = -15)
end
axλ.plot([],[], color = "k", label = "numeric")
axλ.plot([],[], color = "k", label = "model", ls = "dashed")
axλ.legend(ncol = 2)
axgc.set_yticklabels(["", "0.4", "", "", "1.0"])
axλ.set_ylabel("\$\\lambda\$")
axgc.set_ylabel("\$\\kappa / \\kappa(0)\$", labelpad = -20)

# %% Plot a and its uncertainty
axa.clear()
add_grid!(axa);
for i in eachindex(stats)
    Δts, as, δΓs = stats[i]
    p = mean.([log.(a) for a in as])
    axa.plot(Bs2, p, coolcolors[i], label = "\$r=$(rs[ridx[i]])\$")
    # axa.plot(Bs[Bidxs], std.(as), coolcolors[4], label = "\$\\sigma_a\$")
end
axa.legend()
axa.set_xlabel("\$B\$", labelpad = -15)
axa.set_xticks([0, 1])
axa.set_yticks((1,3))
axa.set_ylabel("\$ \\langle \\log(a)\\rangle \$", labelpad = 0)

# %% Plot linear growth and parabolic fits
Δts, as, δΓs = stats[2]
axΓ.clear()
zoomin.clear()
zbox = (0, 0.5)
d = (4, 2)

for (j, ax) in enumerate((axΓ,))
    line, = ax.plot([zbox[1], zbox[1] + d[1], zbox[1] + d[1], zbox[1], zbox[1]],
         [zbox[2], zbox[2], zbox[2]+d[2], zbox[2]+d[2], zbox[2]],
         color = "black", lw = 2.0, zorder = 99)

         line.set_clip_on(false)
end

maxt = 20.0
for (j, i) in enumerate(Bplot)
    B, tc, inc_c = getindex.((Bs2, Δts, δΓs), i)

    axΓ.scatter(tc[1:10:end], exp.(inc_c)[1:10:end],
             color = coolcolors[j], alpha = 0.2, s = 15.0)
    zoomin.scatter(tc, exp.(inc_c),
             color = coolcolors[j], alpha = 0.02, s = 15.0)

    for (p, ax) in enumerate((axΓ, zoomin))
        T = π/B
        t = 0:0.01:T
        f = z_MPSB.(t, B)
        if p == 2
            ax.plot(t, f, color = "k")
        end
        ax.plot(t, f, color = coolcolors[j], ls = "dashed")

    end

    axΓ.scatter([], [], color = coolcolors[j], label = "\$B=$B\$")
end

con = matplotlib.patches.ConnectionPatch(
xyA=zbox .+ (0,d[2]/2) , xyB=zbox .+ (d[1], d[2]/2), coordsA="data", coordsB="data",
axesA=zoomin, axesB=axΓ, color="black", lw = 1.0, zorder = 99
)
zoomin.add_artist(con)

axΓ.set_ylabel("\$|\\delta\\Gamma_{j+1}| \\;/\\; |\\delta\\Gamma_{j}'|\$");
axΓ.set_xlabel("\$\\Delta t_j\$", labelpad = -10);
zoomin.set_xlabel("\$\\Delta t_j\$", labelpad = -10);
axΓ.set_xlim(-0.5, maxt)
axΓ.set_ylim(-0.5, maxt/2)
axΓ.set_yticks([0, maxt/2])
axΓ.set_xticks([0, maxt])
axΓ.axhline(1, ls = "solid", color = "k", lw = 1.0)
zoomin.axhline(1, ls = "solid", color = "k", lw = 1.0)
axΓ.legend(scatterpoints = 1, handletextpad = 0.1, markerscale  = 2,ncol=2)
zoomin.set_ylim(zbox[2], zbox[2] + d[2])
zoomin.set_xlim(zbox[1], zbox[1] + d[1])
zoomin.set_yticks(sort([zoomin.get_ylim()..., 1.0]))
zoomin.set_xticks(zoomin.get_xlim())



# %% Identifiers, spacings, etc
tight_layout()
subplots_adjust(wspace=0.25, hspace = 0.25, left = 0.1,right = 0.97, bottom = 0.06)
add_identifiers!()
# savefig(papersdir()*"Figures/mpsb.png", transparent = true)
