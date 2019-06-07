using DrWatson
include("constants.jl")
quickactivate(@__DIR__, "MagneticBilliardsLyapunovs")
using DynamicalBilliards, Statistics, BSON, LsqFit
using Measurements: measurement, value, uncertainty

include(srcdir()*"mpsb_cellstatistics.jl")

BSON.@load datadir()*"periodic_sinai/Λ.bson" Λ Λσ rs Bs
BSON.@load datadir()*"periodic_sinai/κ.bson" κ κσ

ridx = (3, 5)
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
            tc, inc_c, jhs = manystatistics(r, B, 1000; T = 500)
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
axλ = subplot(3,2,1)
axgc = subplot(3,2,2)
axdist = subplot(3,2,4)
axa = subplot(3,2,3)
axΓ = subplot(3,2,5)
zoomin = subplot(3,2,6)

# %% Toy model
# c1 is always mean collision time
parabola1(x, p) = @. 1.0  + p[1]*x + p[2]*x^2
parabola2(x, p) = @. sqrt(p[1]*sin(x) + p[2]*x + p[3]*x^2)
cosine1(x,B) = @. sqrt(1 + ((1 - cos(2B*x))/2B)^2 + (sin(2B*x)/(2B))^2)
toymodel1(c1, c2, c3) = (log(c2) + log(parabola1(c1, c3)))/c1
toymodel2(c1, c2, c3) = (log(c2) + log(parabola2(c1, c3)))/c1
toymodel3(c1,c2,c3,c4) = (log(c2) - 0.5(c4/c2)^2 + log(parabola1(c1, c3)))/c1
toymodel4(c1,c2,c3) = (c2 + log(parabola1(c1, c3)))/c1
toymodel5(c1, c2, c3) = (c2 + log(cosine1(c1,c3)))/c1

axλ.clear()
params_plot = []
for i in eachindex(ridx)
    Δts, as, δΓs = stats[i]
    params1 = []
    params2 = []

    # for j in eachindex(Bs2)
    #     fit = curve_fit(parabola1, Δts[j], exp.(δΓs[j]), ones(2))
    #     push!(params1, fit.param)
    #     # fit = curve_fit(parabola2, Δts[i], exp.(δΓs[i]), ones(3))
    #     # push!(params2, fit.param)
    # end

    aμ = mean.(as)
    # aσ = std.(as)
    logam = mean.([log.(z) for z in as])
    # toy1 = toymodel1.(κ[Bidxs, ridx[i]], aμ, params1)
    # toy2 = toymodel2.(κ[Bidxs, ridx[i]], aμ, params2)
    # toy3 = toymodel4.(κ[Bidxs, ridx[i]], logam, params1)
    toy4 = toymodel5.(κ[Bidxs, ridx[i]], logam, Bs[Bidxs])

    axλ.plot(Bs, Λ[:, ridx[i]], color = coolcolors[i], label = "r=$(rs[ridx[i]])")
    # coolfill(Bs, Λ[:, ridx[i]], Λσ[:, ridx[i]], axλ, coolcolors[i], "r=$(rs[ridx[i]])")


    # axλ.plot(Bs2, toy1, color = coolcolors[i], ls = "dotted")
    # axλ.plot(Bs2, toy2, color = coolcolors[i], ls = "solid")
    # axλ.plot(Bs2, toy3, color = coolcolors[i], ls = "dashed")
    axλ.plot(Bs2, toy4, color = coolcolors[i], ls = "dotted")
    if i == 1
        push!(params_plot, params1, params2)
    end
end
# axλ.plot([], [], color = "k", ls = "solid", label = "numeric")
# axλ.plot([], [], color = "k", ls = "dashed", label = "toy, \$\\log(\\langle a \\rangle)\$")
# axλ.plot([], [], color = "k", ls = "dotted", label = "toy, \$\\langle \\log(a) \\rangle\$")
axλ.legend(ncol = 2)
axλ.set_ylim(0.5, 2.5)
axλ.set_xlabel("\$B\$", labelpad = -15)
axλ.set_xticks([0, 1])
axλ.set_ylabel("\$\\lambda\$")

# %% Plot a and its uncertainty
axa.clear()
for i in eachindex(stats)
    Δts, as, δΓs = stats[i]
    p = mean.([log.(a) for a in as])
    axa.plot(Bs2, p, coolcolors[i], label = "\$r=$(rs[ridx[i]])\$")
    # axa.plot(Bs[Bidxs], std.(as), coolcolors[4], label = "\$\\sigma_a\$")
end
axa.legend()
axa.set_xlabel("\$B\$", labelpad = -15)
axa.set_xticks([0, 1])
axa.set_yticks(1:2)
axa.set_ylabel("\$ \\langle \\log(a)\\rangle \$", labelpad = 0 )

# %% Histogram of values of a
Δts, as, δΓs = stats[1]
axdist.clear()
axdist_x = 20.0
for (j, i) ∈ collect(enumerate(Bplot))[1:2:end]
    B, jh = getindex.((Bs[Bidxs], as), i)
    coolhist(axdist, jh, 0:0.1:axdist_x, coolcolors[j], "\$B=$B\$")
end
axdist.legend()
axdist.axvline(1, color = "k", lw = 1)
axdist.set_xlim(0, axdist_x)
axdist.set_ylim(0, 0.5)
axdist.set_yticks([])
axdist.set_xlabel("\$a\$", labelpad = -15)
axdist.set_xticks([0, axdist_x])
axdist.set_ylabel("\$P(a)\$", labelpad = 5)


# %% Plot linear growth and parabolic fits
axΓ.clear()
zoomin.clear()


zbox = (0, 0.5)
d = (4, 2)

for (j, ax) in enumerate((axΓ,))
    line, = ax.plot([zbox[1], zbox[1] + d[1], zbox[1] + d[1], zbox[1], zbox[1]],
         [zbox[2], zbox[2], zbox[2]+d[2], zbox[2]+d[2], zbox[2]],
         color = "black", lw = 1.0, zorder = 99)

         line.set_clip_on(false)
end

params1, params2 = params_plot
maxt = 20.0
for (j, i) in enumerate(Bplot)
    B, tc, inc_c = getindex.((Bs2, Δts, δΓs), i)

    axΓ.scatter(tc[1:10:end], exp.(inc_c)[1:10:end],
             color = coolcolors[j], alpha = 0.2, s = 15.0)
    zoomin.scatter(tc, exp.(inc_c),
             color = coolcolors[j], alpha = 0.02, s = 15.0)

    axΓ.scatter([], [], color = coolcolors[j], label = "\$B=$B\$")
    # Parabolic fit:
    tt = 0.0:0.01:maxt
    axΓ.plot(tt, parabola1.(tt, Ref(params1[i])), color = coolcolors[j], ls ="dashed")
    zoomin.plot(tt, parabola1.(tt, Ref(params1[i])), color = coolcolors[j], ls ="dashed")
    # axΓ.plot(tt, parabola2.(tt, Ref(params2[i])), color = coolcolors[j], ls ="dotted")
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
axΓ.legend(scatterpoints = 1, markerscale  = 2, fontsize = 22,ncol=2)
zoomin.set_ylim(zbox[2], zbox[2] + d[2])
zoomin.set_xlim(zbox[1], zbox[1] + d[1])
zoomin.set_yticks(zoomin.get_ylim())
zoomin.set_xticks(zoomin.get_xlim())
# axΓ.set_yticks(0:10:30)
# axΓ.set_xticks(0:10:30)


# %% Identifiers, spacings, etc
tight_layout()
subplots_adjust(wspace=0.2, hspace = 0.2, left = 0.1,right = 0.97, bottom = 0.06)
add_identifiers!()
savefig(papersdir()*"figures/mpsb.png")
