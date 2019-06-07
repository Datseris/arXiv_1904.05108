using DynamicalBilliards


"""
    mushroom_unitcells(obst)
Given the obstacles of the result of `perturbationgrowth` compute the
unitcells of the mushroom billiard.

These are returned as a vector of 3-tuples. First index is start of the cell,
second index is end of the cell and third index is phase (`1` for chaotic,
`2` for laminar).
"""
function mushroom_unitcells(obst)
    L = length(obst)
    ucells = Vector{Tuple{Int, Int, Int}}()
    capobst = (3, 4, 5); cap = 4
    stemobst = (1, 2, 6)

    last_cap = findfirst(x -> x == cap, obst) # first cap collision
    isnothing(last_cap) && return ucells
    i = last_cap + 2

    while i ≤ length(obst) - 2
        if obst[i] ∈ capobst # we start a laminar cell
            final_cap = findfirst(x -> obst[x] ∈ stemobst, i:2:L)
            isnothing(final_cap) && return ucells
            final_cap = (i:2:L)[final_cap] # go to actual indices
            final_cap -= 2 # to go to final cap collision
            push!(ucells, (i - 1, final_cap, 2))
            i = final_cap + 2 # set i to next obstacle
        else # if next obstacle is not in capobst, it must be in stemobst
            next_cap = findfirst(x -> obst[x] == cap, i:2:L)
            isnothing(next_cap) && return ucells
            next_cap = (i:2:L)[next_cap] # go to actual indices
            push!(ucells, (i - 1, next_cap, 1))
            i = next_cap + 2 # set i to next obstacle
        end
    end
    return ucells
end

"""
    psb_unitcells(obst)
Given the obstacles of the result of `perturbationgrowth` compute the
unitcells of the PSB/MPSB billiard.

These are returned as a vector of 2-tuples. First index is start of the cell,
second index is end of the cell.
"""
function psb_unitcells(obst)
    L = length(obst)
    ucells = Vector{Tuple{Int, Int}}()
    disk = 1

    isempty(obst) && return ucells

    i = findfirst(x -> obst[x] == disk, obst) # first disk collision
    isnothing(i) && return ucells
    i += 2 # first obstacle after disk

    while i ≤ length(obst) - 2 # increment 2, because all obst indices are doubled
        next_disk = findfirst(x -> obst[x] == disk, i:2:L)
        isnothing(next_disk) && return ucells
        next_disk = (i:2:L)[next_disk] # go to actual indices
        push!(ucells, (i - 1, next_disk))
        i = next_disk + 2 # set i to next obstacle
    end
    return ucells
end
