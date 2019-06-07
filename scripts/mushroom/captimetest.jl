using DrWatson
quickactivate(@__DIR__, "MagneticBilliardsLyapunovs")
using BSON, PyPlot
l = 1.0
data = BSON.load(datadir()*"mushrooms/mushroom_capviastem_l=$(l)_r=1.0_T=10000.0_N=10000.bson")

plot(data[:ws], data[:t_c])
w = data[:ws]

plot(w, 2 .* (@. 1 - w^2/(36) - w^4/(1200) - w^6/(15680) - w^8/(145152)) .+ π*l)

#=
import DynamicalBilliards:MushroomTools
using SymPy

Y(r, w) = 1/2π * (sqrt(4r^2 - w^2) * SymPy.SpecialFuncs.elliptic_e( w^2/(-4r^2 + w^2) ) + π*r * Float64( SymPy.mpmath[:hyp3f2]( 0.5, 0.5, 0.5, 1, 1.5, w^2/(4r^2) ) ) )
a(r, w)  = SymPy.mpmath[:hyp3f2](0.5, 0.5, 0.5, 1, 1.5, w^2/(4r^2))

YY = 2 .* Y.(1.0, w) .+ π*l

plot(w, YY)
=#
