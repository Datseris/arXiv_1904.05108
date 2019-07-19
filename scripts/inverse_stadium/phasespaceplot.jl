using Revise, DrWatson
quickactivate(@__DIR__, "MagneticBilliardsLyapunovs")
using InverseStadium, DynamicalBilliards, Makie, StatsMakie

to_θ(vx, vy) =  @. mod(atan(vy/vx) + (vx < 0)*π, 2π)

function gen_data(ω = 10.0, N = 100, T = 1000.0,dt = 0.01)
    bds = billiard_dual_stadium()
    xs = Float64[]
    ys = Float64[]
    θs = Float64[]

    for i ∈ 1:N
        ps = dual_inside(bds, ω)
        dual_bounce!(ps, bds)

        x,y, vx, vy = timeseries(ps, bds, T, dt)

        θ =  to_θ(vx, vy)

        append!(xs, x)
        append!(ys, y)
        append!(θs, θ)
    end
    return xs, ys, θs
end

function full3Dplot(xs, ys, θs)
    scn = Scene()
    # careful: memory is finite!
    scatter!(scn, xs,ys,θs, markersize = 0.01, color = θs, colormap=:phase, overdraw = true)
    return scn
end

function section2D(xs, ys, θs, f)
    toplot = Point2{Float64}[]
    colors = Float64[]
    for (x,y,θ) in zip(xs, ys, θs)
        plotpoint, xp, yp, cp =  f(x,y,θ)
        if plotpoint
            push!(toplot, Point2(xp, yp))
            push!(colors, cp)
        end
    end

    scn = Scene()

    scatter!(scn, toplot, color = colors, colormap=:phase, markersize = 0.01,  alpha = 0.1)
    return scn
end

function hist2D(xs, ys, θs, f)
    plotx = Float64[]
    ploty = Float64[]
    colors = Float64[]
    for (x,y,θ) in zip(xs, ys, θs)
        plotpoint, xp, yp, cp =  f(x,y,θ)
        if plotpoint
            push!(plotx, xp)
            push!(ploty, yp)
            push!(colors, cp)
        end
    end

    scn = Scene()

    heatmap!(scn, histogram(nbins = (50,100)), plotx, ploty, colormap=:viridis)
    return scn
end


function psos_θ(θ_c = 0.0; δ = 0.01)
    function f(x,y,θ)
        if abs(θ - θ_c) < δ
            return true, x, y, θ
        else
            return false, 0.0, 0.0, 0.0
        end
    end
    return f
end

xs, ys, θs = gen_data(2.0, 500, 100.0, 0.001)

hist2D(xs, ys, θs, psos_θ(0, δ = 2π))
