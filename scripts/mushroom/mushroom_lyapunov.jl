using DrWatson
quickactivate(@__DIR__, "MagneticBilliardsLyapunovs")
using DynamicalBilliards, Statistics

const N = 5000
ws = 0.05:0.05:2
hs = 0.05:0.05:2

Λ = zeros(length(ws), length(hs))
Λσ = zeros(length(ws), length(hs))

for (j, w) in enumerate(ws)
    println("w = $w"); flush(stdout)
    @time for (i, h) in enumerate(hs)
        bd = billiard_mushroom(h, w, 1.0)
        ps = [MushroomTools.randomchaotic(h,w,1.0) for _ in 1:N]
        # λs = parallelize(lyapunovspectrum, bd, N, ps)
        λs = [lyapunovspectrum!(p, bd, N)[1] for p in ps]
        Λ[i, j] = mean(λs)
        Λσ[i, j] = std(λs)
    end
end

description =
"""
Λ is the mean value of the exponent.
Λσ is the standard deviation.
The first indices are the `hs`, the second are the `ws`.
"""

using BSON

BSON.bson(datadir()*"mushrooms/Λ_N=$N.bson", (@dict Λ Λσ ws hs description))
