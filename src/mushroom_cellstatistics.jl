include("unitcells.jl")
using LinearAlgebra

function statistics_tc_tl(ts, Δ, obst)
    ucells = mushroom_unitcells(obst)

    tc = Float64[]; tl = Float64[]; tlc = Float64[]
    inc_c = Float64[]; inc_l = Float64[]; inc_lc = Float64[]
    jhs = Float64[]

    for i in 1:length(ucells)-1
        exp(Δ[ucells[i][1]]) == Inf || exp(Δ[ucells[i][2]]) == Inf && break
        s, f = ucells[i]
        δt = ts[f] - ts[s]
        if i > 1 && ucells[i-1][3] == 2 && ucells[i][3] == 1 # chaotic after laminar
            push!(tlc, δt)
            push!(inc_lc, Δ[f] - Δ[s])
        elseif ucells[i][3] == 1 # chaotic cell
            push!(tc, δt)
            push!(inc_c, Δ[f] - Δ[s])
            if !isnan(exp(Δ[f+1] - Δ[f])) && !isinf(exp(Δ[f+1] - Δ[f]))
                push!(jhs, exp(Δ[f+1] - Δ[f]))
            end
        else
            push!(tl, δt)
            push!(inc_l, Δ[f] - Δ[s])
        end
    end
    if length(tl) ∉ (length(tlc), length(tlc) + 1)
        println(length(tl), " ", length(tlc))
        error("For every laminar one chaotic")
    end
    return tc, inc_c, tl, inc_l, tlc, inc_lc, jhs
end

function manystatistics_tc_tl(l, w, N = 2000, T = 5000.0)
    bd = billiard_mushroom(l,w,1.0)
    ttc = Float64[]; ttl = Float64[]
    tinc_c = Float64[]; tinc_l = Float64[]
    ttlc = Float64[]; tinc_lc = Float64[]
    tjs = Float64[]

    for n in 1:N
        p = MushroomTools.randomchaotic(l,w,1.0)
        ts, Gs, obst = perturbationgrowth(p, bd, T)
        Δs = log.(norm.(perturbationevolution(Gs)))
        tc, inc_c, tl, inc_l, tlc, inc_lc, jh = statistics_tc_tl(ts, Δs, obst)
        for (a, x) in zip((tc, inc_c, tl, inc_l, tlc, inc_lc, jh), (ttc, tinc_c, ttl, tinc_l, ttlc, tinc_lc, tjs))
            append!(x, a)
        end
    end
    return (ttc, tinc_c, ttl, tinc_l, ttlc, tinc_lc, tjs)
end;
