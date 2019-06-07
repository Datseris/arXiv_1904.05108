using StaticArrays
using LinearAlgebra

import DynamicalBilliards:lyapunovspectrum, lyapunovspectrum!
export lyapunovspectrum, lyapunovspectrum!


"""
    boundary_transition
applies the matrix T from the Vörös paper to the 2D offset vectors
"""
function boundary_transition!(offset::Vector{SVector{2, T}}, p::AbstractParticle{T},
                             Δω::T, o::Obstacle) where {T}

    sφ = DynamicalBilliards.cross2D(p.vel, - DynamicalBilliards.normalvec(o, p.pos))
    tanφ = sφ/sqrt(1-sφ^2) # tan(arcsin(x)) = x/sqrt(1 - x^2)

    for k in 1:2
        δΓ = offset[k]
        offset[k] = SVector{2,T}(δΓ[1], δΓ[2] + Δω*tanφ*δΓ[1])
    end
end


function propagate_offset2D!(offset::Vector{SVector{2, T}}, t::T,
                            p::Particle{T}) where {T}
    for k in 1:2
        δΓ = offset[k]
        offset[k] = SVector{2,T}(δΓ[1] + t*δΓ[2], δΓ[2])
    end

end


function propagate_offset2D!(offset::Vector{SVector{2, T}}, t::T,
                            p::MagneticParticle{T}) where {T}
    s, c = sincos(p.omega*t)

    for k in 1:2
        δΓ = offset[k]
        offset[k] = c*δΓ + s*SVector{2,T}(δΓ[2]/p.omega, -p.omega*δΓ[1])
    end

end

"""
    lyapunovspectrum(ps::DualParticle, bds::DualBilliard, t)
Lyapunov spectrum for inverse stadium billiards, uses the coordinates described by
Vörös et al. because they have already calculated all the tangent maps…

"""
lyapunovspectrum(ps::DualParticle, args...) = lyapunovspectrum!(deepcopy(ps), args...)



function lyapunovspectrum!(ps::DualParticle{T}, bds::DualBilliard{T},
                                tt::AbstractFloat) where {T}

    offset = [SVector{2, T}(1,0), SVector{2, T}(0,1)]

    t = T(tt)

    count = zero(T)
    λ = zeros(T, 2)

    ω = ps[2].omega

    while count < t
        #bounce linear particle
        i, tmin, cp = next_collision(ps[1], bds[1])

        o = bds[1][i]
        relocate_outside!(ps[1], o, tmin, cp)
        propagate_offset2D!(offset, tmin, ps[1])
        
        #handle first boundary transition
        boundary_transition!(offset, ps[1], ω, o)
        
        update!(ps[2], ps[1])
        
        #increment counter
        count += DynamicalBilliards.increment_counter(t, tmin)
        
        #bounce magnetic particle
        i, tmin, cp = next_collision(ps[2], bds[2])

        if tmin == T(Inf)
            throw(ErrorException("pinned particle"))
        end

        o = bds[2][i]
        relocate_outside!(ps[2], o, tmin, cp)
        ps[2].center = find_cyclotron(ps[2])

        propagate_offset2D!(offset, tmin, ps[2])

        boundary_transition!(offset, ps[2], -ω, bds[2][i])

        update!(ps[1], ps[2])

        count += DynamicalBilliards.increment_counter(t, tmin)

        # decompose & profit
        
        Q, R = qr(hcat(offset[1], offset[2]))
        offset[1], offset[2] = Q[:, 1], Q[:, 2]
        for i ∈ 1:2
            λ[i] += log(abs(R[i,i]))
        end

    end

    return λ./count
end


## experimental function, please ignore ##############################################
#=
#individual component matrices
P(t::T) where T = SMatrix{2,2, T}( 1,  t,
                                   0  ,1 )

E(t::T, ω) where T = SMatrix{2,2, T}(cos(ω*t),     1/ω * sin(ω*t),
                                     -ω*sin(ω*t),  cos(ω*t)      )

Ta(Δω::T, φ) where T =  SMatrix{2,2,T}(1,          0,
                                       Δω*tan(φ),  1)

#product
It(ω, t_lin, t_mag, φ_in, φ_out) = Ta(-ω, φ_in)*E(t_mag, ω)*Ta(ω, φ_out)*P(t_lin)

"""
    get_stability
Produces modern art
"""
function get_stability(bds,ω, Nξ::Int, Nsφ::Int = Nξ)
    δξ = totallength(bds[1])/Nξ
    δφ = 2/Nsφ

    d = Dict{SV{Int}, Float64}()


    # I don't have another easy way to create a dual particle
    # I'll have to add one…
    dummy = dual_inside(bds, ω)


    for ξcell ∈ 1:Nξ, φcell ∈ 1:Nsφ
        #get center position
        ξc = (ξcell - 0.5)*δξ
        φc = (φcell - 0.5)*δφ - 1


        #convert to real space
        pos, vel, i = from_bcoords(ξc, φc, bds[1])

        #set dummy coordinates
        dummy[1].pos = pos
        dummy[1].vel = vel
        update!(dummy[2], dummy[1])

        τ, s1, s2 = dual_bounce!(dummy, bds)

        φ_in = asin(DynamicalBilliards.cross2D(s2[4], -normalvec(bds[2][s2[1]], s2[3])))
        φ_out = asin(DynamicalBilliards.cross2D(s1[4], -normalvec(bds[1][s1[1]], s1[3])))

        DEBUG && println("ω = $(ω)\n t_lin = $(s1[2])\n t_mag = $(s2[2])\n φ_in = $(φ_in)\n φ_out = $(φ_out)")

        it = It(ω, s1[2], s2[2], φ_in, φ_out)

        d[SV{Int}(ξcell,φcell)] = tr(it)
    end
    return d
end
=#
