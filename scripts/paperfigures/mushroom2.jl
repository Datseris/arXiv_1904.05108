using DrWatson
include("constants.jl")
quickactivate(@__DIR__, "MagneticBilliardsLyapunovs")
using DynamicalBilliards, Statistics, BSON, LsqFit
using Measurements: measurement, value, uncertainty

include(srcdir()*"mushroom_cellstatistics.jl")

import FileIO

f = BSON.load(datadir()*"/mushrooms/Λ_N=5000.bson")
Λ = f[:Λ]; Λws = f[:ws]; Λls = f[:ls]


file = datadir()*"/mushrooms/unit_cell_times.bson"
BSON.@load file result hs ws
hi = 10; hhh = hs[hi];
tc, tl = result[hi, :, 1], result[hi, :, 2]
tlc, fl = result[hi, :, 3], result[hi, :, 4]
καππα(h, w) = MushroomTools.g_c_3D(h, w, 1.0) * (π*(h*w + π/2))/w


# # %%
# wval = 1.6; lval = 1.0
# tc, inc_c, tl, inc_l, tlc, inc_lc, jh = manystatistics_tc_tl(lval, wval, 500);
# figure()
# scatter(tc, exp.(inc_c));
# ylabel("\$|\\delta\\Gamma|_n \\;/\\; |\\delta\\Gamma|_{n-1}\$")
# scatter(tlc, exp.(inc_lc));
# scatter(tl, exp.(inc_l)); xlabel("\$t\$");
# title("w = $wval"); tight_layout();

# %% Compute parameters of toy model
using LsqFit: curve_fit
linearfit(x, p) = @. p[1] + p[2]*x
lin(τ, a, b) = a + τ*b;

function toyparameters(l, w, N = 2000, T = 2000.0)
    tc, inc_c, tl, inc_l, tlc, inc_lc, jh = manystatistics_tc_tl(l, w, N, T)

    logj = mean(log.(jh))
    as = Float64[]; bs = Float64[]
    τs = Float64[]

    lf = length(tl)/(length(tc)+length(tl)) # laminar phase frequency

    for (t, inc) in zip((tc, tl, tlc), (inc_c, inc_l, inc_lc))
        xx = findall(isinf, inc)
        deleteat!(t, xx)
        deleteat!(inc, xx)
        xx = findall(isnan, inc)
        deleteat!(t, xx)
        deleteat!(inc, xx)
        τ = mean(t)
        push!(τs, τ)
        fit = curve_fit(linearfit, t, exp.(inc), [1.0, 1.0])
        a, b = fit.param
        push!(as, a)
        push!(bs, b)
    end

    logs = [mean(log(as[i] + bs[i]*t) for t in (tc, tl, tlc)[i] if as[i] + bs[i]*t > 0) for i in 1:3]

    # # Linear fit of combination of laminarchaotic
    # tll = [tl[i] + tlc[i] for i in 1:min(length(tl), length(tlc))]
    # inc_ll = [inc_l[i] + inc_lc[i] for i in 1:min(length(inc_l), length(inc_lc))]
    # τ = mean(tll)
    # push!(τs, τ)
    # fit = curve_fit(linearfit, tll, exp.(inc_ll), [1.0, 1.0])
    # a, b = fit.param
    # push!(as, a)
    # push!(bs, b)

    # Linear fit of all possible chaotic phases
    tcc = vcat(tc, tlc)
    inc_cc = vcat(inc_c, inc_lc)
    τcc = mean(tcc)
    push!(τs, τcc)
    fit = curve_fit(linearfit, tcc, exp.(inc_cc), [1.0, 1.0])
    a, b = fit.param
    push!(as, a)
    push!(bs, b)

    # Treat laminar and chaotic after laminar as ONE episode
    inc_2 = inc_l[1:length(inc_lc)] .* inc_lc
    t2 = tl[1:length(tlc)] .+ tlc

    push!(logs, mean(log(as[4] + bs[4]*t) for t in  tcc if as[4] + bs[4]*t > 0))
    push!(logs, mean(log(as[4] + bs[4]*t ) for t in  tcc if as[4] + bs[4]*t > 0))

    return logj, lf, τs, as, bs, logs
end

# %% Gather parameters versus w
# collect parameters
WTEST = HTEST = 0.1:0.1:2.0
HS = WS = [0.5, 1.0, 1.5]
N = 2000; T = 5000.0
filename = datadir()*"mushrooms/toypars_N=$(N)_T=$(T).bson"

if !isfile(filename)
    toypar_w = [[] for l in HS]
    for (z, h) in enumerate(HS)
        println("h = $h")
        for (i, w) in enumerate(WTEST)
            toyp = toyparameters(h, w, N, T)
            push!(toypar_w[z], toyp)
        end
    end

    # Similarly collect versus h
    toypar_h = [[] for l in HS]
    for (wi, w) in enumerate(WS)
        println("w = $w")
        for h in HTEST
            toyp = toyparameters(h, w, N, T)
            push!(toypar_h[wi], toyp)
        end
    end
    BSON.@save filename toypar_h toypar_w
else
    fff = BSON.load(filename)
    toypar_h = fff[:toypar_h]
    toypar_w = fff[:toypar_w]
end


# %%
figure(figsize = (2figx, 10))
axλ = subplot(2,2,1)
axλh = subplot(2,2,3)
axΓ = subplot(2,4,3)
zoomin = subplot(2,4,7)
axaw = subplot(2,4,4)
axah = subplot(2,4,8)

# figure()
# axλ = subplot(211)
# axλh = subplot(212)


# %% Lyapunovs versus w:
axλ.clear()
axλh.clear()

function toy(logj, lf, τs, as, bs, logs, κ)
    nom = logj + log(
                   (1-lf)  *  lin(τs[1], as[1], bs[1]) +
                   (lf) * lin(τs[2], as[2], bs[2]) +
                   (lf) * lin(τs[3], as[3], bs[3])
                )
    T = lf*(τs[2] + τs[3]) + (1-lf)*τs[1]
    T = κ
    return nom/T
end
function toymerged(logj, lf, τs, as, bs, logs, κ)
    nom = logj + log(
                    lin(τs[4], as[4], bs[4]) +
                   (lf) * lin(τs[2], as[2], bs[2])
                )
    # T = lf*(τs[2]) + τs[1]
    T = κ
    return nom/T
end

function logavrg(logj, lf, τs, as, bs, logs, κ)
    nom = logj + (1-lf)*logs[1] + lf*(logs[2] + logs[3])
    T = κ
    return nom/T
end
function mergedc(logj, lf, τs, as, bs, logs, κ)
    nom = logj + logs[4] + lf*(logs[2])
    T = κ
    # T = lf*(τs[2]) + τs[1]
    return nom/T
end

toymodels = (toymerged,)
toplot = [1,3]






# actual plotting here:

for hi in toplot
    h = HS[hi]

    lx = findfirst(isequal(h), Λls)
    λs = Λ[lx, :]
    axλ.plot(Λws, λs, color = coolcolors[hi], label = "\$h = $h\$")


    toyλ = [zeros(length(WTEST)) for _ in toymodels]
    for wi in eachindex(toypar_w[hi])
        toypar = toypar_w[hi][wi]
        for (fi, f) in enumerate(toymodels)
            toyλ[fi][wi] = f(toypar..., καππα(h, WTEST[wi]))
        end
    end
    for (fi, f) in enumerate(toymodels)
        if fi == 1
            axλ.plot(HTEST, toyλ[fi], ls = "dashed",
            color = coolcolors[hi]);
        else
            axλ.plot(HTEST, toyλ[fi], ls = "dotted",
            color = coolcolors[hi], marker = markers[fi]);
        end
    end
end
axλ.legend(ncol = 2, fontsize = 22)

for wi in toplot
    w = WS[wi]

    lx = findfirst(isequal(w), Λls)
    λs = Λ[:, lx]
    axλh.plot(Λws, λs, color = coolcolors[wi], label = "\$w = $w\$")


    toyλ = [zeros(length(HTEST)) for _ in toymodels]
    for hi in eachindex(toypar_h[wi])
        toypar = toypar_h[wi][hi]
        for (fi, f) in enumerate(toymodels)
            toyλ[fi][hi] = f(toypar..., καππα(HTEST[hi], w))
        end
    end
    for (fi, f) in enumerate(toymodels)
        if fi == 1
            axλh.plot(HTEST, toyλ[fi], ls = "dashed",
            color = coolcolors[wi]);
        else
            axλh.plot(HTEST, toyλ[fi], ls = "dotted",
            color = coolcolors[wi], marker = markers[fi]);
        end
    end
end

# Legends:
axλh.legend(ncol = 2, fontsize = 22)
axλh.set_ylabel("\$\\lambda\$", labelpad = -15)
axλ.set_ylabel("\$\\lambda\$", labelpad = -15)

for ax in (axλ, axλh); ax.set_xticks([0, 2]); end
axλ.set_yticks([0.3, 0.5])
axλh.set_yticks([0.2, 0.4])
axλ.set_xlabel("\$w\$", labelpad = -20)
axλh.set_xlabel("\$h\$", labelpad = -20)

# %% Plot toymodel parameters and theoretical curves
idx = 2 # which value of w or h
@assert WTEST == HTEST
@assert HS == WS

axah.clear()
axaw.clear()
for (ax, toypar, sym) in zip((axaw, axah), (toypar_w[idx], toypar_h[idx]), (:w, :h))

    # Plot other parameters:
    loga = [x[1] for x in toypar]
    ττs = [x[3] for x in toypar]

    # more accurate τs:
    if sym == :w
        τs = result[findfirst(isequal(HS[idx]), hs), :, 1:3]
        fl = result[findfirst(isequal(HS[idx]), hs), :, 4]
    else
        τs = result[:, findfirst(isequal(WS[idx]), ws), 1:3]
        fl = result[:, findfirst(isequal(WS[idx]), ws), 4]
    end

    for (i, c) in enumerate(["c", "\\ell"])
        if i == 2
            ax.plot(WTEST, τs[:, i]/10, #[τ[i]/10 for τ in τs],
            label = "\$\\tau_{$(c)}/10\$", color = coolcolors[i+3])
        else
            ax.plot(WTEST, (1 .- fl) .* τs[:, 1] + fl .* τs[:, 3], #
            label = "\$\\tau_{$(c)}\$", color = coolcolors[i+3])
            # ax.plot(WTEST, [t[4] for t in ττs])
        end
        # ax.plot(WTEST, [a[i] for a in as], color = coolcolors[i+3],
        #      ls = "dotted", lw = 1.0, marker = "D")
        # ax.plot(WTEST, [b[i] for b in bs], color = coolcolors[i+3],
        #      ls = "dotted", lw = 1.0, marker = "o")
    end
    ax.plot(WTEST, 10fl, color = coolcolors[6], label="\$10f_\\ell\$")
    ax.plot(WTEST, 5loga, color = coolcolors[7], label="\$5\\langle\\log(a)\\rangle\$")

    ax.set_xlabel("\$$sym\$")
    # ax.axhline(0, color = "k", lw = 1.0)
end

for ax in (axaw, axah); ax.set_xticks([0, 2]); end
axaw.set_xlabel("\$w\$", labelpad = -20)
axah.set_xlabel("\$h\$", labelpad = -20)
axaw.set_ylabel("parameters", labelpad = 0)
axah.set_ylabel("parameters", labelpad = 0)
leg = axah.legend(ncol = 2, fontsize = 22)

# set the linewidth of each legend object
for legobj in leg.legendHandles
    legobj.set_linewidth(4.0)
end

# %% Linear fit in δΓ vs t
axΓ.clear()
zoomin.clear()
zbox = (0, -1)
d = (10, 10)

for (j, ax) in enumerate((axΓ,))
    line, = ax.plot([zbox[1], zbox[1] + d[1], zbox[1] + d[1], zbox[1], zbox[1]],
         [zbox[2], zbox[2], zbox[2]+d[2], zbox[2]+d[2], zbox[2]],
         color = "black", lw = 2.0, zorder = 99)

         line.set_clip_on(false)
end

con = matplotlib.patches.ConnectionPatch(
xyA=zbox .+ (d[1]/2,d[2]) , xyB=zbox .+ (d[1]/2, 0), coordsA="data", coordsB="data",
axesA=zoomin, axesB=axΓ, color="black", lw = 1.0, zorder = 99
)
zoomin.add_artist(con)


z_MB(t) = sqrt(1 + 0.9^2*t^2 - 1/2*t)

let hi = 3, wi = 15 # wi = length(WTEST)÷2
    ttt = 0:0.01:100.0
    h = HS[hi]
    w = WTEST[wi]
    tc, inc_c, tl, inc_l, tlc, inc_lc, jh = manystatistics_tc_tl(h, w, 500)
    for (j, (t, inc)) in enumerate(zip((tc, tl, tlc), (inc_c, inc_l, inc_lc)))
        for f in (isinf, isnan)
            xx = findall(f, inc)
            deleteat!(t, xx)
            deleteat!(inc, xx)
        end
        la = ["c", "\\ell", ""][j]

        for (fff, ax) in enumerate((axΓ, zoomin))
            ax.scatter(t, exp.(inc), color = coolcolors[j < 3 ? j+3 : j+1],
            alpha = fff == 1 ? 0.25 : 0.05, s = fff == 1 ? 20 : 15)
            logj, lf, τs, as, bs,  = toypar_w[hi][wi]
            a = as[j]; b = bs[j]

            if j < 3
                ax.scatter([], [], color = coolcolors[j < 3 ? j+3 : j+1],
                s = 5, label = "\$$(la)\$")
            end

            if j == 2
                ax.plot(ttt, lin.(ttt, a, b),
                color = coolcolors[j+3], ls = "dashed")
            end
            # ax.plot(ttt, z_MB.(ttt),color = coolcolors[6], ls = "dotted")
        end

        # ax.plot(ttt, z_MB.(ttt),color = coolcolors[6], ls = "dotted")
    end
    for (fff, ax) in enumerate((axΓ, zoomin))
        logj, lf, τs, as, bs,  = toypar_w[hi][wi]
        ax.plot(ttt, lin.(ttt, as[4], bs[4]),
        color ="black", ls = "dashed")
    end

    axΓ.text(0.1, 0.8, "\$h=$h,w=$w\$", transform = axΓ.transAxes,
    va = "center")
end

axΓ.set_ylim((0,50) .+ zbox[2])
axΓ.set_xlim(0,50)
axΓ.set_xticks(axΓ.get_xlim())
axΓ.set_yticks(axΓ.get_ylim())
axΓ.set_ylabel("\$|\\delta\\Gamma_{j+1}| \\;/\\; |\\delta\\Gamma_{j}'|\$", labelpad = -20);
# axΓ.set_xlabel("\$\\Delta t_j\$", labelpad = -20);

zoomin.set_ylabel("\$|\\delta\\Gamma_{j+1}| \\;/\\; |\\delta\\Gamma_{j}'|\$", labelpad = -20);
zoomin.set_ylim(zbox[2], zbox[2] + d[2])
zoomin.set_xlim(zbox[1], zbox[1] + d[1])
zoomin.set_yticks(sort([zoomin.get_ylim()..., 0.0]))
zoomin.set_xticks(sort([zoomin.get_xlim()...]))
zoomin.set_xlabel("\$\\Delta t_j\$", labelpad = -20);
zoomin.legend(scatterpoints = 1, markerscale  = 4,fontsize = 22,ncol=1,
framealpha=1.0,handletextpad=0.1)
zoomin.axhline(0, color = "k", lw = 1.0)

# %%
tight_layout()
add_identifiers!()
subplots_adjust(wspace=0.25, hspace = 0.2, left = 0.04,right = 0.97, bottom = 0.08)

savefig(papersdir()*"figures/mushroom2.png", transparent = true)
