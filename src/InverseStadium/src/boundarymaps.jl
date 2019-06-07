import DynamicalBilliards:boundarymap, boundarymap_portion, phasespace_portion, SV
export boundarymap, boundarymap_portion, phasespace_portion

"""
    boundarymap(bds, t, ps)
Boundary map for dual particles. Returns
* a tuple containing the arclengths and angles of incidence for particles
  coming from the inside
* the arclength intervals of the first billiard
"""
function boundarymap(pss::AbstractVector{<:DP{T}}, bds::DB{T}, t,
                     intervals = arcintervals(bds[1])) where {T}

    bmaps = Vector{SV{T}}[]
    for p ∈ pss
        bm, i = boundarymap(p, bds, t, intervals)
        push!(bmaps, bm)
    end

    return bmaps, intervals
end
        


function boundarymap(pars::DP{T}, bds::DB{T}, t,
                     intervals = arcintervals(bds[1])) where {T}

    ps = deepcopy(pars)
    bd = bds[1]
    bmap = SV{T}[]
    if typeof(t) == Int
        sizehint!(bmap, t)
    end
    count = zero(t); t_to_write = zero(T)

    while count < t
        tmin, state1 = dual_bounce!(ps,bds)
        t_to_write += tmin
        i = state1[1]
        
        ξ, sφ = to_bcoords(state1[3], state1[4], bd[i])
 
        @inbounds push!(bmap, SV{T}(ξ + intervals[i], sφ))
        # set counter
        count += DynamicalBilliards.increment_counter(t, t_to_write)
        t_to_write = zero(T)
    end #time, or collision number, loop
    return bmap, intervals
end


################################################################################
## Boundary map & phase space portions
################################################################################

function boundarymap_portion(bds::DB{T}, t, pars::DP{T}, δξ, δφ = δξ;
                             intervals = arcintervals(bds[1])) where {T}
    ps = deepcopy(pars)
    count = zero(T)
    t_to_write = zero(T)

    d = Dict{SV{Int}, Int}()
    bd = bds[1]
    p = ps[1]

    while count < t
        tmin, state1, = dual_bounce!(ps,bds)
        t_to_write += tmin

        i = state1[1]
        if typeof(bd[i]) <: PeriodicWall
            continue # do not write output if collision with with PeriodicWall
        else
            # get Birkhoff coordinates
            ξ, sφ = to_bcoords(state1[3], state1[4], bd[i])
            ξ += intervals[i]
            # compute index & increment dictionary entry
            ind = SV{Int}(floor(Int, ξ/δξ), floor(Int, (sφ + 1)/δφ))

            d[ind] = get(d, ind, 0) + 1

            # set counter
            count += DynamicalBilliards.increment_counter(t, t_to_write)
            t_to_write = zero(T)
        end
    end #time or collision number loop

    #calculate ratio of visited boxes
    total_boxes = ceil(Int, totallength(bd)/δξ) * ceil(Int, 2/δφ)
    ratio = length(keys(d))/total_boxes

    return ratio, d
end

function boundarymap_portion(bds::DB{T}, t, parss::AbstractVector{<:DP{T}}, δξ, δφ = δξ;
                             intervals = arcintervals(bds[1])) where {T}    
    d = Dict{SV{Int}, Int}()
    bd = bds[1]

    intervals = arcintervals(bd)
    for pars in parss
        count = zero(T)
        t_to_write = zero(T)

        ps = deepcopy(pars)
        p = ps[1]

        while count < t
            tmin, state1, = dual_bounce!(ps,bds)
            t_to_write += tmin

            i = state1[1]
            if typeof(bd[i]) <: PeriodicWall
                continue # do not write output if collision with with PeriodicWall
            else
                # get Birkhoff coordinates
                ξ, sφ = to_bcoords(state1[3], state1[4], bd[i])
                ξ += intervals[i]
                # compute index & increment dictionary entry
                ind = SV{Int}(floor(Int, ξ/δξ), floor(Int, (sφ + 1)/δφ))

                d[ind] = get(d, ind, 0) + 1

                # set counter
            end
                count += DynamicalBilliards.increment_counter(t, t_to_write)
                t_to_write = zero(T)
        end #time or collision number loop
    end

    #calculate ratio of visited boxes
    total_boxes = ceil(Int, totallength(bd)/δξ) * ceil(Int, 2/δφ)
    ratio = length(keys(d))/total_boxes

    return ratio, d
end


function phasespace_portion(bds::DB{T}, pars::DP{T}, dict, δξ, δφ = δξ) where {T}
    bd = bds[1]
    
    ints = arcintervals(bd)

    maxξ = ceil(Int, totallength(bd)/δξ)
    maxφ = ceil(Int, 2/δφ)

    dummy = deepcopy(pars)

    total = zero(T); visited = zero(T)

    
    for ξcell ∈ 1:maxξ-1, φcell ∈ 1:maxφ-1
        #get center position
        ξc = (ξcell - 0.5)*δξ
        φc = (φcell - 0.5)*δφ - 1

        pos, vel, i = from_bcoords(ξc, φc, bd, ints)

        dummy[1].pos = pos
        dummy[1].vel = vel
        DynamicalBilliards.specular!(dummy[1], bd[i]) # because vel has to point outwards!
        update!(dummy[2], dummy[1])
        τ,s1,s2 = dual_bounce!(reverse(dummy),reverse(bds))

        total += τ
        (haskey(dict, SV{Int64}(ξcell, φcell))) && (visited += τ)

    end

    return visited/total, visited, total#, dict
end


phasespace_portion(bds::DB, t, pars::DP, δξ, δφ = δξ; intervals = arcinterval(bds[1])) = 
    phasespace_portion(bds, pars, boundarymap_portion(bds, t, pars, δξ, δφ, intervals = intervals)[2], δξ, δφ)
  
