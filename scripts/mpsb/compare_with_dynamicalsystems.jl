using DynamicalSystems, DrWatson, PyPlot, OrdinaryDiffEq
quickactivate(@__DIR__, "MagneticBilliardsLyapunovs")

diffeq = (alg = Vern9(), dtmax = 0.05, reltol = 1e-9, abstol = 1e-9,
maxiters = 1e9)

B = 1.0; c = 0.1; d0 = 0.5
r = d0/2

ds = Systems.antidots(rand(4); B = B, c = c, d0 = d0)
p = ds.p

function antidot_init(p; x = rand(), y = rand(), θ = -π + 2*π*rand())
    var0 = zeros(Float64, 4)
    var0[1] = x
    var0[2] = y
    # Check if the initial condition is valid:
    pot0 = Systems.antidot_potential(var0[1], var0[2], p)
    while pot0 >= 0.99
        var0[1] = rand()
        var0[2] = rand()
        pot0 = Systems.antidot_potential(var0[1], var0[2], p)
    end
    θ0 = θ
    s = sqrt(2*(1 - pot0))
    var0[3] = s*cos(θ0)
    var0[4] = s*sin(θ0)
    return var0
end

# %%
figure()
u0 = antidot_init(p)
tr = trajectory(ds, 100.0, u0; diffeq...)
x, y, vx = columns(tr)
limits = [minimum(x), maximum(x), minimum(y), maximum(y)]

#point density for potential plot
n = 512
xlin = range(limits[1], limits[2], length = n)
ylin = range(limits[3], limits[4], length = n)
# define arrays necesarry for plotting
xgrid = repeat(xlin', n, 1)
ygrid = repeat(ylin ,1,n)
pot = zeros(Float64, n, n)

# calculate potential for plotting
for col in 1:n, row in 1:n
    pot[row, col] = Systems.antidot_potential(xlin[col], ylin[row], p)
end
ax = gca()
cpfilled = ax[:contourf](xgrid, ygrid, pot, 8, alpha=1.0, vmin = 0.1,
extend="max", cmap="Greys", levels=0:0.25:1.25)
# add red contour line at potential = 1.0
c1line = ax[:contour](xgrid, ygrid, pot, colors="red", linewidths=1.0,
levels = [1.0])
plot(x, y)













# %% Compare lyapunovs without B
using BSON
BSON.@load datadir()*"periodic_sinai/Λ0.bson" Λ0 rs

ds = Systems.antidots(rand(4); B = 0.0, c = c, d0 = d0)

x0, y0, θ0 = 0.5, 0.5, π/4

nl = 1
u0 = antidot_init(p; x = x0, y = y0, θ = θ0)
tinteg = tangent_integrator(ds, nl; diffeq...)
# pinteg = parallel_integrator(ds, [u0, u0 + 1e-9rand(SVector{4})]; diffeq...)
λs = []
for r in rs
    set_parameter!(ds, 2, 2r)
    λ = zeros(nl); X = 0
    @time for i in 1:10
        ic = SVector{4}(antidot_init(ds.p))
        reinit!(tinteg, ic, orthonormal(4, nl))
        l = lyapunovs(tinteg, 5000, 1.0, 10.0)
        if l[1] > 0.1
            λ .+= l
            X += 1
        else
            println("not large")
        end
    end
    push!(λs, λ./X)
    # reinit!(pinteg)
    # @time push!(λs, lyapunov(pinteg, 2000.0, 10.0, 0.1, 1e-9, 1e-6, 1e-12))
end

figure()
plot(rs, Λ0, label = "billiard")
plot(rs, λs, label = "ode")
xlabel("r"); ylabel("λ"); legend()










# %% compare for magnetic fields
using BSON
r = 0.15
BSON.@load datadir()*"periodic_sinai/Λ.bson" Λ rs Bs
ri = findfirst(isequal(r), rs)

bb = 0:0.01:1.0

ds = Systems.antidots(rand(4); B = 0.0, c = c, d0 = 2*rs[ri])

nl = 1
tinteg = tangent_integrator(ds, nl; diffeq...)
# pinteg = parallel_integrator(ds, [u0, u0 + 1e-9rand(SVector{4})]; diffeq...)
λs = []
for B in bb
    println(B)
    set_parameter!(ds, 1, B)
    λ = zeros(nl); n = 0
    @time for i in 1:10
        ic = SVector{4}(antidot_init(ds.p))
        reinit!(tinteg, ic, orthonormal(4, nl))
        l = lyapunovs(tinteg, 5000, 1.0, 10.0)
        if l[1] > 0.1
            λ .+= l
            n += 1
        else
            println("not large")
        end
    end
    push!(λs, λ./n)
    # reinit!(pinteg)
    # @time push!(λs, lyapunov(pinteg, 2000.0, 10.0, 0.1, 1e-9, 1e-6, 1e-12))
end

figure()
subplot(211)
plot(Bs, Λ[:, ri], label = "billiard")
plot(bb, λs ./ sqrt(2), label = "ode")
xlabel("B"); ylabel("λ"); legend()
title("r = $r")
subplot(212)
plot(Bs, (λs[2:end] ./ sqrt(2)) .- Λ[:, ri], label = "difference")
xlabel("B"); ylabel("λ"); legend()
title("r = $r")

using BSON
BSON.bson(datadir()*"periodic_sinai/Λ_antidots_r=$(r).bson", @dict bb λs)
