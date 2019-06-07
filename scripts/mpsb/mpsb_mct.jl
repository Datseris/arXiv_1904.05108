using DrWatson
quickactivate(@__DIR__, "MagneticBilliardsLyapunovs")
using DynamicalBilliards, Statistics

const N = 5000
Bs = 0.01:0.05:1
rs = 0.1:0.05:0.45

κ = zeros(length(Bs), length(rs))
κσ = zeros(length(Bs), length(rs))

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
        # λs = parallelize(meancollisiontime, bd, N, ps)
        κs = [meancollisiontime!(p, bd, N) for p in ps]
        κ[j, i] = mean(κs)
        κσ[j, i] = std(κs)
    end
end

description =
"""
κ is the mean value of the m.c.t. for different B and r in the MPSB.
κσ is the standard deviation.
The first indices are the `Bs`, the second are the `rs`.
"""
using BSON
BSON.bson(datadir()*"periodic_sinai/κ_N=$N.bson", (@dict κ κσ Bs rs description))

# %% Same process but without B
κ0 = zeros(length(rs)); κ0σ = copy(κ0)

for (i, r) in enumerate(rs)
    bd = billiard_sinai(r; setting = "periodic")
    println("r = $r")
    ps = [randominside(bd) for _ in 1:N]
    # λs = parallelize(lyapunovspectrum, bd, N, ps)
    λs = [meancollisiontime(p, bd, N) for p in ps]
    κ0[i] = mean(λs)
    κ0σ[i] = std(λs)
end

BSON.bson(datadir()*"periodic_sinai/κ0_N=$N.bson", (@dict κ0 κ0σ rs))
