"obstacle color based on index"
function mushcolor(index::Int)
    if index ∈ (2,6) # stem sides
        return "C0"
    elseif index == 1 # stem bottom
        return "C2"
    elseif index == 4 # cap semicircle
        return "C3"
    elseif index ∈ (3,5) # cap bottom
        return "C1"
    else
        return "k"
    end
end
sinaicolor(x) = x == 1 ? "C3" : "C0"

"obstacle point size based on index"
mushsize(x) = x ∈ (2,3,5,6) ? 0.6 : 1.2
sinaisize(x) = x == 1 ? 1.0 : 0.5

"""
    colorcodeplot(ts, Δ, obst)
"""
function colorcodeplot(ts, logΔ, obst;
    tlim = (0, maximum(ts)), ax = gca(),
    obstsize = mushsize, obstcolor = mushcolor)

    i, e = t_to_i.(tlim, Ref(ts))
    plotrange = i:e

    # adaptive sizing of points
    relsize = 20sqrt(length(ts)/(e - i))

    ax.scatter(ts[plotrange], logΔ[plotrange], color = obstcolor.(obst[plotrange]),
            s = relsize .* obstsize.(obst[plotrange]), zorder = 10)
    ax.plot(ts[plotrange], logΔ[plotrange], color="k", lw=0.5, zorder = 0)
    return plotrange
end

"convert collision time `t` to index in `ts`"
function t_to_i(val, A)
    i = 1
    d = abs(val - A[i])
    @inbounds for j in eachindex(A)
        dd = abs(val - A[j])
        if dd < d
            i = j
            d = dd
        end
    end
    return i
end
