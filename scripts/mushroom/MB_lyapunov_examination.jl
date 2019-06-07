using DrWatson
quickactivate(@__DIR__, "MagneticBilliardsLyapunovs")
using DynamicalBilliards, Statistics, BSON, LsqFit, Parameters
using Measurements: measurement, value, uncertainty

include(srcdir()*"mushroom_cellstatistics.jl")

import FileIO

f = BSON.load(datadir()*"/mushrooms/Λ_N=5000.bson")
Λ = f[:Λ]; Λws = f[:ws]; Λls = f[:ls]

include(scriptdir()*"paperfigures/constants.jl")
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
linear01fit(x, p) = @. 1.0 + p[1]*x
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

    # Linear fit of all possible chaotic phases
    tcc = vcat(tc, tlc)
    inc_cc = vcat(inc_c, inc_lc)
    push!(τs, mean(tcc)) # 4th time
    fit = curve_fit(linearfit, tcc, exp.(inc_cc), [1.0, 1.0])
    a, b = fit.param
    push!(as, a)
    push!(bs, b)
    push!(logs, mean(log(as[4] + bs[4]*t) for t in  tcc if as[4] + bs[4]*t > 0))

    # Assume all chaotic episodes have laminar, even of length zero
    M = length(tc) - length(tl)
    @assert M > 0
    append!(tl, zeros(M))
    append!(inc_l, zeros(M))
    fit = curve_fit(linear01fit, tl, exp.(inc_l), [1.0])
    b, = fit.param
    push!(as, 1.0)
    push!(bs, b)
    push!(τs, mean(tl))
    push!(logs, mean(log(as[5] + bs[5]*t) for t in tl if as[5] + bs[5]*t ≥ 0))

    return logj, lf, τs, as, bs, logs
end

# %% Gather parameters
WTEST = HTEST = 0.1:0.1:2.0
HS = WS = [0.5, 1.0, 1.5]

function g(d)
    @unpack N, T = d

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

    return @dict toypar_w toypar_h
end

N = 2000; T = 2000.0
file = produce_or_load(datadir()*"mushrooms/toytest", @dict(N, T), g; force = true)
@unpack toypar_w, toypar_h = file

# %%
using PyPlot
figure(figsize = (18, 8))
axλ = subplot(1,2,1)
axλh= subplot(1,2,2)


# %% Lyapunovs versus w:
axλ.clear()
axλh.clear()

function toymerged(logj, lf, τs, as, bs, logs, κ)
    nom = logj + log(
                    lin(τs[4], as[4], bs[4]) +
                   (lf) * lin(τs[2], as[2], bs[2])
                )
    # T = lf*(τs[2]) + τs[1]
    T = κ
    return nom/T
end

function alwayslaminar(logj, lf, τs, as, bs, logs, κ)
    nom = logj + log(
                    lin(τs[4], as[4], bs[4]) +
                    lin(τs[5], as[5], bs[5])
                )
    # T = lf*(τs[2]) + τs[1]
    T = κ
    println("a = ", as[5])
    println("b = ", bs[5])
    println("τ = ", τs[5])
    println("τ2 = ", τs[2]*lf)

    return nom/T
end


toymodels = (toymerged,alwayslaminar)
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
