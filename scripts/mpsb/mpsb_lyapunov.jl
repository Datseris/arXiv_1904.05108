using DrWatson
quickactivate(@__DIR__, "MagneticBilliardsLyapunovs")
using DynamicalBilliards, Statistics

const N = 5000
Bs = 0.01:0.005:1
rs = 0.1:0.05:0.45

Λ = zeros(length(Bs), length(rs))
Λσ = zeros(length(Bs), length(rs))

function randomparticles(bd, B)
    ps = MagneticParticle{Float64}[]
    while length(ps) < N
        p = randominside(bd, 2B)
        !ispinned(p, bd) && push!(ps, p)
    end
    return ps
end

for (i, r) in enumerate(rs)
    bd = billiard_sinai(r; setting = "periodic")
    println("r = $r"); flush(stdout)
    @time for (j, B) in enumerate(Bs)
        ps = randomparticles(bd, B)
        # λs = parallelize(lyapunovspectrum, bd, N, ps)
        λs = [lyapunovspectrum!(p, bd, N)[1] for p in ps]
        Λ[j, i] = mean(λs)
        Λσ[j, i] = std(λs)
    end
end

description =
"""
Λ is the mean value of the exponent for different B and r in the MPSB.
Λσ is the standard deviation.
The first indices are the `Bs`, the second are the `rs`.
"""

using BSON

BSON.bson(datadir()*"periodic_sinai/Λ_N=$N.bson", (@dict Λ Λσ Bs rs description))

# %% Same process but without B
Λ0 = zeros(length(rs)); Λ0σ = copy(Λ0)

for (i, r) in enumerate(rs)
    bd = billiard_sinai(r; setting = "periodic")
    println("r = $r")
    ps = [randominside(bd) for _ in 1:N]
    # λs = parallelize(lyapunovspectrum, bd, N, ps)
    λs = [lyapunovspectrum(p, bd, N)[1] for p in ps]
    Λ0[i] = mean(λs)
    Λ0σ[i] = std(λs)
end

BSON.bson(datadir()*"periodic_sinai/Λ0_N=$N.bson", (@dict Λ0 Λ0σ rs))
