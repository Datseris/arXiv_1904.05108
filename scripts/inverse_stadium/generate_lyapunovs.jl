using Revise, DrWatson
quickactivate(@__DIR__, "MagneticBilliardsLyapunovs")
using InverseStadium, DynamicalBilliards, PyPlot, Printf

function single_particle_lyapunov(;l = 1.0, w = 1.0, ω = 1.0, T = 1000.0)
    bds = billiard_dual_stadium(l, w)
    ps = dual_inside(bds, ω)
    λ = lyapunovspectrum(ps, bds, Float64(T))
    return ps, λ[1]
end


function many_particle_lyapunovs(N = 1000, T = 1000.0; l = 1.0, w = 1.0, ω = 1.0, save = false)
    ps = Vector{DP}(undef, N)
    λs = Vector{Float64}(undef, N)

    pinned = 0

    for i ∈ 1:N
        got_result = false
        while !got_result
            try
                ps[i], λs[i] = single_particle_lyapunov(l = l, w = w, ω = ω, T = T)
                got_result = true
            catch e
                if !isa(e, InterruptException)
                    pinned += 1
                else
                    throw(e)
                end
            end
        end
    end

    if save == true
        bson("data/isb_lyap_lω=$(@sprintf("%.03f", log10(ω)))_N=$(N)_T=$(T)_l=$(l)_w=$(w).bson",
             Dict(:particles => ps, :lyapunovs => λs, :pinned => pinned, :ω => ω))
    end
    return ps, λs, pinned
end


function multithreaded_lyapunovs(ωs;N = 1e4, T = 1e4, l = 1.0, w = 1.0)

    Threads.@threads for ω ∈ ωs
        many_particle_lyapunovs(N, T, l = l, w = w, ω = ω, save = true)
    end

end

#multithreaded_lyapunovs(10 .^(0:0.125:3), N = 1000, T = 1e4)
