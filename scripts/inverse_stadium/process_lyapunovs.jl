using Revise, DrWatson
quickactivate(@__DIR__, "MagneticBilliardsLyapunovs")
using InverseStadium, DynamicalBilliards, PyPlot, BSON, Statistics, StaticArrays

"""
    filter_regular(pss, λss)
filters out "regular" particles by cutting all Lyapunov exponents below an 
empirical threshold, which is `0.1*maximum(λss)`
"""
function filter_regular(pss, λss)
    ps::Vector{DP{Float64}} = copy(pss)
    λs = copy(λss)

    threshold = 1/10 * maximum(λs) # empirical value

    del_ind = Int[]
    for i ∈ 1:length(λs)
        if λs[i] < threshold
            push!(del_ind, i)
        end
    end

    # delete regular particles
    deleteat!(ps, del_ind)
    deleteat!(λs, del_ind)

    return ps, λs, length(del_ind)
end

"""
    get_psr(ps, l=1.0, w=1.0; T = 10000.0, Nboxes = 200)
Returns a value proportional (?) to the phase space volume taken up by `ps` and
a value proportional to the total phase space volume.
"""
function get_psr(ps, l=1.0, w=1.0; T = 10000.0, Nboxes = 200)
    bds = billiard_dual_stadium(l, w)

    δξ = totallength(bds[1])/Nboxes
    δφ = 2/Nboxes

    _, dict = boundarymap_portion(bds, T, ps, δξ, δφ)

    # needs one particle as template for dummy, but uses data from all particles
    _, psr, pst = phasespace_portion(bds, ps[1], dict, δξ, δφ)
    
    
    return psr, pst
end

"""
    eval_data()
calls all the functions on all the data in the `data`-subdirectory.
"""
function eval_data(Tbox  = 100.0, Nbox = 100, save = true)
    files = readdir("data/")

    ωs = Vector{Float64}(undef, length(files))
    λss = copy(ωs)
    σλs = copy(ωs)
    psv_vis = copy(ωs)
    psv_tot = copy(ωs)
    println("Tbox = $(Tbox); Nbox = $(Nbox)")
    Threads.@threads for i in 1:length(files)

        println("Thread $(Threads.threadid()): "*files[i])
        dict = BSON.load("data/"*files[i])       
        ps, λs = filter_regular(dict[:particles], dict[:lyapunovs])
        # technically, I would have to add arguments for l and w here, but we always use 1.0 SO FAR!
        visited, total = get_psr(ps, T = Tbox, Nboxes = Nbox)
        
        ωs[i] = dict[:ω]
        λss[i] = mean(λs)
        σλs[i] = std(λs)
        psv_vis[i] = visited
        psv_tot[i] = total
    end
    save && bson("output_Tbox=$(Tbox)_Nbox=$(Nbox).bson", @dict ωs λss σλs psv_vis psv_tot)
    return ωs, λss, σλs, psv_vis, psv_tot
end

ωs, λss, σλs, psv_vis, psv_tot = eval_data(10000.0, 200)
