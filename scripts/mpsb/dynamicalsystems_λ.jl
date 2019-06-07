using DrWatson
quickactivate(@__DIR__, "MagneticBilliardsLyapunovs")
using DynamicalSystems, OrdinaryDiffEq

diffeq = (alg = Vern9(), dtmax = 0.05, reltol = 1e-9, abstol = 1e-9,
maxiters = 1e9)

# %% compare for magnetic fields
using BSON
BSON.@load datadir()*"periodic_sinai/Λ.bson" Λ rs Bs

for r in (0.15, 0.3)
    if isfile(datadir()*"periodic_sinai/Λ_antidots_r=$(r).bson")
        continue
    end
    ds = Systems.antidots(rand(4); B = 0.01, c = 0.1, d0 = 2*rs[ri])

    nl = 1
    tinteg = tangent_integrator(ds, nl; diffeq...)
    # pinteg = parallel_integrator(ds, [u0, u0 + 1e-9rand(SVector{4})]; diffeq...)
    λs = Float64[]
    for B in Bs
        println(B)
        set_parameter!(ds, 1, B)
        λ = zeros(nl); n = 0
        @time for i in 1:10
            ic = SVector{4}(antidot_init(ds.p))
            reinit!(tinteg, ic, orthonormal(4, nl))
            l = lyapunovs(tinteg, 5000, 1.0, 10.0)
            if l[1] > 0.1
                λ .+= l
                n += 1
            else
                println("not large")
            end
        end
        push!(λs, λ[1]/n)
    end

    BSON.bson(datadir()*"periodic_sinai/Λ_antidots_r=$(r).bson", @dict Bs λs)
end

# %% plot differences
using PyPlot
figure()
for r in (0.15, 0.3)
    BSON.@load datadir()*"periodic_sinai/Λ_antidots_r=$(r).bson" Bs λs
    ri = findfirst(isequal(r), rs)
    d = abs.(((λs ./ sqrt(2)) .- Λ[:, ri]))
    plot(Bs, d, label = "r = $r")
end
legend();
title("abs diff. of Billiard and DynamicalSystems.jl")
xlabel("B"); ylabel("λ diff."); tight_layout()
