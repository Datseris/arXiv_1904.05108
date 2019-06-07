using DrWatson
quickactivate(@__DIR__, "MagneticBilliardsLyapunovs")
include("constants.jl")
using DynamicalBilliards, Statistics, BSON
using Measurements: measurement, value, uncertainty

include(srcdir()*"mpsb_cellstatistics.jl")

BSON.@load datadir()*"periodic_sinai/Λ0.bson" Λ0 Λ0σ rs
BSON.@load datadir()*"periodic_sinai/κ0_N=2000.bson" κ0 κ0σ
ridx = 1:2:8

# Gather statistics
if !isdefined(Main, :Δts) || length(Δts) < 1
    Δts = Vector{Float64}[];
    as = copy(Δts); δΓs = copy(Δts)
    for (i, r) in enumerate(rs)
        println(i/length(rs))
        tc, inc_c, jhs = manystatistics(r, 0, 5000; T = 2000)
        push!(Δts, tc); push!(as, jhs); push!(δΓs, inc_c)
    end
end

# %% Make figure
figure(figsize = (12, 12))
axλ = subplot(2,2,1)
axa = subplot(2,2,2)
axlinear = subplot(2,2,3)
axdist = subplot(2,2,4)

# %% Toy model
z_PSB(t) = sqrt(1 + (t)^2)
# c1 is always mean collision time
# most basic sqrt + 1 model:
toymodel1(c1, c2, c3) = (log(c2) + log(sqrt(1 + c1^2)))/c1
# averages of logarithms
# toymodel2(c1, c2, c3) = (c2 + c3)/c1
toymodel3(c1, c2, c3) = (c2 + log(z_PSB(c1)))/c1

aμ = mean.(as)
aσ = std.(as)
am = measurement.(aμ, aσ)
κm = measurement.(κ0, κ0σ)
toym = toymodel1.(κm, am, 5)

using Elliptic
b = @. sqrt((0.5 + κ0^2/2)/(1.0 + κ0^2/2)) # dq/dp = 1
b = @. sqrt((κ0^2)/(1 + κ0^2)) # dq/dp = 0
# b = @. sqrt((0.0 + κ0^2/2)/(0.5 + κ0^2/2))
A = @. b*π/rs
theory = @. sqrt(1 + A^2)
theory2 = @. sqrt(1 + A^2 - 2A)

theory3 = @. sqrt(4*b^2 + rs^2)*acsch(2b/rs)/rs + log(b/rs)

logas = [log.(z) for z in as]
logaμ = mean.(logas)
logaσ = std.(logas)

axλ.clear()
# coolfill(rs, Λ0, Λ0σ, axλ, coolcolors[1], "\$\\lambda\$: numeric")
axλ.plot(rs, Λ0, color = coolcolors[1], label = "numeric")

axλ.plot(rs, toymodel3.(κ0, logaμ, 1), color = coolcolors[2], ls = "dashed",
label = "model, numeric")

axλ.plot(rs, toymodel3.(κ0, theory3, 1), color = coolcolors[3],
label = "model, analytic")

# axλ.plot(rs, toymodel1.(κ0, aμ, 1), color = coolcolors[2],
# label = "toy, \$\\log \\left( \\langle a \\rangle \\right)\$",
# ls = "dotted")

# axλ.plot(rs, toymodel1.(κ0, theory, 1), color = coolcolors[4],
# label = "theory 1", ls = "dotted")

axλ.legend()
axλ.set_yticks([1, 2, 3, 4])
axλ.set_xlabel("\$r\$")
axλ.set_ylabel("\$\\lambda\$")

# %% Plot a and its uncertainty versus theoretical value
axa.clear()

# axa.plot(rs, aμ, coolcolors[3], label = "\$\\langle a\\rangle\$, numeric")
# # axa.plot(rs, aσ, coolcolors[4], label = "\$\\sigma_a\$")
# axa.plot(rs, theory, coolcolors[5], label = "\$\\langle a\\rangle\$, analytic")
# axa.plot(rs, theory2, coolcolors[6], ls = "dashed", label = "theory 2")
# axa.plot(rs, exp.(theory3), label = "theory 3")


axa.plot(rs, logaμ, coolcolors[2], ls = "dashed",
label = " numeric")
axa.plot(rs, theory3, coolcolors[3], ls = "solid",
label = "analytic")
axa.legend()
axa.set_yticks(1:3)
axa.set_xlabel("\$r\$")
axa.set_ylabel("\$\\langle \\log(a)\\rangle\$")




# %% Plot linear growth
axlinear.clear()

zbox = (0, 0)
d = 4
# Make zoom-in plot
zoomin = axdist
zoomin.clear()

for (j, ax) in enumerate((axlinear,))
    line, = ax.plot([zbox[1], zbox[1] + d, zbox[1] + d, zbox[1], zbox[1]],
         [zbox[2], zbox[2], zbox[2]+d, zbox[2]+d, zbox[2]],
         color = "black", lw = 1.0, zorder = 99)

         line.set_clip_on(false)
end
# for j in (0, 1)
con = matplotlib.patches.ConnectionPatch(
xyA=zbox .+ (0,d/2) , xyB=zbox .+ (d, d/2), coordsA="data", coordsB="data",
axesA=zoomin, axesB=axlinear, color="black", lw = 1.0, zorder = 99
)
zoomin.add_artist(con)
# end
# zoomin.axis("off")
zoomin.set_xticks([zbox[1], d])
# zoomin.set_yticks(zbox[2]:(zbox[2]+d))
zoomin.set_yticks([zbox[2], 1.0, d])
zoomin.set_xlim(zbox[1], zbox[1] + d)
zoomin.set_ylim(zbox[2], zbox[2] + d)

# plot linear approximation
for ax in (zoomin,)
    ax.plot([0, 100], [0, 100],
        lw = 1.0, color = "k"
    )
    # plot square root approximation
    (z = 0:0.01:100)
    ax.plot(z, z_PSB.(z), ls = "dashed",
    color = "k", lw = 3.0)
end

for i in ridx
    r, tc, inc_c = getindex.((rs, Δts, δΓs), i)
    axlinear.scatter(tc[1:10:end], exp.(inc_c)[1:10:end],
             color = coolcolors[i÷2 + 1], alpha = 0.2, s = 15.0)
    zoomin.scatter(tc[1:10:end], exp.(inc_c)[1:10:end],
             color = coolcolors[i÷2 + 1], alpha = 0.2, s = 10.0)
    axlinear.scatter([], [], color = coolcolors[i÷2 + 1], label = "\$r=$r\$")
end

axlinear.set_ylabel("\$|\\delta\\Gamma_{j+1}| \\;/\\; |\\delta\\Gamma_{j}'|\$");
axlinear.set_xlabel("\$\\Delta t_j\$");
zoomin.set_xlabel("\$\\Delta t_j\$");
axlinear.set_xlim(-1, 25)
axlinear.set_ylim(-1, 25)
axlinear.legend(scatterpoints = 1, markerscale  = 2,  handletextpad = 0.1)
axlinear.set_yticks(0:10:30)
axlinear.set_xticks(0:10:30)

# %% Identifiers, spacings, etc
tight_layout()
subplots_adjust(wspace=0.25, hspace = 0.25, left = 0.1,right = 0.97, bottom = 0.1)
add_identifiers!()
# savefig(papersdir()*"figures/psb.png", transparent = true)
