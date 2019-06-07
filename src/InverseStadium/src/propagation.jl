import DynamicalBilliards:SV, increment_counter, find_cyclotron, totallength,
    timeseries, timeseries!

export dual_inside, dual_evolve, dual_evolve!, dual_bounce!, timeseries, timeseries!


"""
    update!(w::AbstractParticle, r::AbstractParticle)
Updates the state of `w` to match that of `r`
"""
update!(w::AbstractParticle, r::AbstractParticle) = _update!(w, r)

function update!(w::MagneticParticle, r::AbstractParticle)
    _update!(w, r)
    w.center = find_cyclotron(w)
end

@inline function _update!(w::AbstractParticle, r::AbstractParticle)
    w.pos = r.pos
    w.vel = r.vel
    w.current_cell = r.current_cell
    return nothing
end


################################################################################
## actual propagation
################################################################################


"""
    relocate_outside!(p, o, tmin)
Like [`relocate!`](@ref) but ensures the particle is outside the billiard
"""
function relocate_outside!(p::AbstractParticle{T}, o::Obstacle{T}, tmin::T, cp::SV{T}) where {T}
    sig = +1
    propagate!(p, cp, tmin) # propagate to collision point
    d = DynamicalBilliards.distance(p.pos, o)

    # ensure particle is outside of billiard regardless of obstacle type
    okay = d ≤ 0.0
    if okay
        n = DynamicalBilliards.normalvec(o, p.pos)
        p.pos -= d*n
    end
    return okay
end

"""
    bounce_outside!(p, bd)
Like [`bounce!`](@ref) but uses [`relocate_outside!`](@ref) instead of [`relocate!`](@ref)
"""
function bounce_outside!(p::AbstractParticle{T}, bd::Billiard{T}) where {T}
    i::Int, tmin::T, cp::SV{T} = next_collision(p, bd)
    if tmin != T(Inf)
        o = bd[i]
        relocate_outside!(p, o, tmin, cp)
    end
    typeof(p) <: MagneticParticle && (p.center = find_cyclotron(p))   
    return i, tmin, p.pos, p.vel
end


## High-levelish functions #####################################################

"""
    dual_bounce!(ps, bds)
Alternately [`bounce_outside!`](@ref)s the particles in `ps` on the billiard tables
in `bds`, [`update!`](@ref)ing the particles in between bounces. Returns
the total colllision time as well as the return values of both individual `bounce`s.
"""
function dual_bounce!(ps::DP{T}, bds::DB{T}; update = true) where {T}
    o1 = bounce_outside!(ps[1], bds[1])

    update!(ps[2], ps[1])

    o2 = bounce_outside!(ps[2], bds[2])
    update && update!(ps[1], ps[2])

    return o1[2]+o2[2], o1, o2
end


"""
    dual_evolve!(ps::Tuple{Particle, MagneticParticle}, bds, t)
[`evolve!`](@ref)s the "dual particles" on the billiard tables using
 [`dual_bounce`](@ref).

Returns collision times, positions, velocities and angular velocities.
"""
function dual_evolve!(ps::DP{T}, bds::DB{T}, t) where {T <: AbstractFloat}
    if t ≤ 0
        throw(ArgumentError("`evolve!()` cannot evolve backwards in time."))
    end

    rt = T[]; push!(rt, 0)
    rpos = SVector{2,T}[]; push!(rpos, ps[1].pos)
    rvel = SVector{2,T}[]; push!(rvel, ps[1].vel)
    rω = T[]; push!(rω, 0)
    count = zero(t)

    while count < t
        tmin, state1, state2 = dual_bounce!(ps, bds)

        push!(rω, 0)          ; push!(rω, ps[2].omega)
        push!(rt, state1[2])  ; push!(rt, state2[2])
        push!(rpos, state1[3]); push!(rpos, state2[3])
        push!(rvel, state1[4]); push!(rvel, state2[4])

        count += increment_counter(t, tmin)
    end

    return rt, rpos, rvel, rω
end


dual_evolve(p::DP, args...) = dual_evolve!(deepcopy(p), args...)


function timeseries!(ps::DP{T}, bds::DB{T}, t, dt = T(0.01),
                     warning::Bool = true) where {T}

    ts = [zero(T)]
    x  = [ps[1].pos[1]]; y  = [ps[1].pos[2]]
    vx = [ps[1].vel[1]]; vy = [ps[1].vel[2]]

    prevω = ps[2].ω
    
    count = zero(t)
    t_total = zero(T)
    t_to_write = zero(T)

    prevpos = ps[1].pos + ps[1].current_cell
    prevvel = ps[1].vel
       
    @inbounds while count < t

        ct, state1, state2 = dual_bounce!(ps, bds, update = false)
        states = (state1, state2)
        
        for p_idx ∈ 1:2 
            i = states[p_idx][1]
            t_to_write += states[p_idx][2]
            
            if DynamicalBilliards.isperiodic(i, bds[p_idx])
                # do nothing at periodic obstacles
                continue
            else
                if t_to_write ≤ dt
                    # push collision point only
                    push!(ts, t_to_write)
                    push!(x, ps[p_idx].pos[1] + ps[p_idx].current_cell[1])
                    push!(y, ps[p_idx].pos[2] + ps[p_idx].current_cell[2])

                    push!(vx, ps[p_idx].vel[1])
                    push!(vy, ps[p_idx].vel[2])
                else
                    # extrapolate & append
                    nx, ny, nvx, nvy, nts =
                        DynamicalBilliards.extrapolate(ps[p_idx], prevpos,
                                                       prevvel, t_to_write,
                                                       dt, prevω)

                    append!(ts, nts[2:end] .+ t_total)
                    append!(x, nx[2:end])
                    append!(y, ny[2:end])
                    append!(vx, nvx[2:end])
                    append!(vy, nvy[2:end])
                end

                prevpos = ps[p_idx].pos + ps[p_idx].current_cell
                prevvel = ps[p_idx].vel
                
                t_total += t_to_write
                count += DynamicalBilliards.increment_counter(t, t_to_write)
                t_to_write = zero(T)

            end
        end
        update!(ps[1], ps[2])
        
    end
    
    return x, y, vx, vy, ts
end

timeseries(p::DP, args...) = timeseries!(deepcopy(p), args...)
