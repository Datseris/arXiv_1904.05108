include("unitcells.jl")
using LinearAlgebra

function statistics(ts, Δ, obst)
    ucells = psb_unitcells(obst)
    tc = Float64[]; inc_c = Float64[]; jhs = Float64[]
    isempty(ucells) && return tc, inc_c, jhs

    for i in 1:length(ucells)-1
        s, f = ucells[i]
        δt = ts[f] - ts[s]
        d = Δ[f] - Δ[s]
        isinf(exp(d)) || isnan(exp(d)) && return tc, inc_c, jhs
        push!(tc, δt)
        push!(inc_c, d)
        if i > 1
            d = exp(Δ[s] - Δ[s-1])
            isinf(d) || isnan(d) && return tc, inc_c, jhs
            push!(jhs, d)
        end
    end
    deleteat!(jhs, findall(isinf, jhs))
    return tc, inc_c, jhs
end

"""
    manystatistics(r, B, N = 1000)
Return vectors of the cell durations, the corresponding logarithmic
increase of perturbation vector during the cell and the jump heights.

Notice that the jump heights are the ratio of the norms before and after.
They don't use logarithms.
"""
function manystatistics(r, B, N = 1000; T = 1000.0)
    bd = billiard_sinai(r; setting = "periodic")
    ttc = Float64[];
    tinc_c = Float64[];
    tjs = Float64[]

    for n in 1:N
        if B == 0
            p = randominside(bd)
        else
            p = randominside(bd, 2B)
        end
        ts, Gs, obst = perturbationgrowth(p, bd, T)
        isempty(ts) && continue
        Δs = log.(norm.(perturbationevolution(Gs)))
        tc, inc_c, jh = statistics(ts, Δs, obst)
        for (a, x) in zip((tc, inc_c, jh), (ttc, tinc_c, tjs))
            append!(x, a)
        end
    end
    return (ttc, tinc_c, tjs)
end;
