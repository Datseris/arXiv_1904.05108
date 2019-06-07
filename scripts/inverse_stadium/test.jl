using Revise, DrWatson
quickactivate(@__DIR__, "MagneticBilliardsLyapunovs")
using InverseStadium
################################################################################
## TEST FUNCTIONS
################################################################################

function dual_test(; ω = 2.0, Nξ = 100, Nφ = Nξ)

    bds = billiard_dual_stadium(1.0, 1.0)
    bd = bds[1]

    ints = arcintervals(bd)

    δξ = totallength(bd)/Nξ; δφ = 2/Nφ


    dummy = dual_inside(bds, ω)

    λs = Matrix{Float64}(undef, Nξ, Nφ)

    pinned = 0
    for ξcell ∈ 1:Nξ, φcell ∈ 1:Nφ
        #get center position
        ξc = (ξcell - 0.5)*δξ
        φc = (φcell - 0.5)*δφ - 1

        pos, vel, i = from_bcoords(ξc, φc, bd, ints)

        dummy[1].pos = pos
        dummy[1].vel = vel
        DynamicalBilliards.specular!(dummy[1], bd[i]) # because vel has to point outwards!
        InverseStadium.update!(dummy[2], dummy[1])
        try
            λs[ξcell, φcell] = InverseStadium.lyapunovspectrum(dummy, bds, 2000.0)[1]
        catch
            pinned += 1
            #println("($(ξc), $(φc)) pinned")
            λs[ξcell, φcell] = NaN
        end

    end

    plt.pcolormesh(λs')
    plt.colorbar()
    println("$(pinned) pinned trajectories")
    return λs, pinned
end

using Statistics, LinearAlgebra

function get_λ_hist(λs= dual_test(); bins = 100)
    cutoff = (0.1 * maximum(λs .* .!(isnan.(λs))))
    println("Cutoff value: $cutoff")
    mask = λs .> cutoff
    λ_masked = λs .* mask
    N = sum(mask)

    µ = sum(λ_masked)/N
    σ = sqrt(1/(N + 1) * sum(@. mask*(λ_masked - µ)^2))

    figure("hist")
    y, x′ = plt.hist((λs.*mask)[:], bins = range(cutoff, maximum(λs.*mask), length=bins))
    x = x′[2:end] .- diff(x′)

    figure()
    plt.plot(x, normalize(y))
    plt.plot(x, normalize(@.(exp(-((x - µ)^2)/(2σ^2)))))
    plt.axvline(µ, color="C2")
end

λs, pinned = dual_test(ω = 10.0)
get_λ_hist(λs)


#λs_ = map(x->(isnan(x)?missing:x), λs)
