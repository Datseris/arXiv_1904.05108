using DrWatson
quickactivate(@__DIR__, "MagneticBilliardsLyapunovs")
using DynamicalBilliards, PyPlot, Statistics
SV = SVector{2}

sloc = 0.0
w = 0.5 # stem width
h = 1.0 # stem length
r = 1.0 # cap radius
w2 = 0.2 # stem width 2
h2 = 0.2
r2 = 0.3

@assert w2 + 2r2 < h
@assert w2 < 2r2

abs(sloc) + w/2 > r && error("Stem is outside the mushroom cap!")

leftcorn = SV(-w/2 + sloc, 0)
rightcorn = SV(w/2 + sloc, 0)
upleftcorn = SV(-w/2 + sloc, h)
uprightcorn = SV(w/2 + sloc, h)

stembot = FiniteWall(leftcorn, rightcorn, SV(0, w), false, "Stem bottom")
stemleft = FiniteWall(upleftcorn, leftcorn, SV(w, 0), false, "Stem left")


farleft = SV(-r, h)
farright = SV(r, h)

capbotleft = FiniteWall(
farleft, upleftcorn, SV(0, w), false, "Cap bottom left")
capbotright = FiniteWall(
uprightcorn, farright, SV(0, w), false, "Cap bottom right")

cap = Semicircle([0.0, h], r, [0.0, -1.0], "Mushroom cap")

# Make space for the other mushroom:
stemrightlow = FiniteWall([w/2, 0], [w/2, h/2 - w2/2], SV(-1, 0), false, "Stem right")
stemrightupper = FiniteWall([w/2, h/2 + w2/2], [w/2, h], SV(-1, 0), false, "Stem right")

# Make the other mushroom:
ob1 = FiniteWall([w/2, h/2-w2/2], [w/2+h2, h/2-w2/2], SV(0, 1))
ob5 = FiniteWall([w/2, h/2+w2/2], [w/2+h2, h/2+w2/2], SV(0, -1))
ob2 = FiniteWall([w/2+h2, h/2-w2/2],[w/2+h2, h/2-r2], SV(1, 0))
ob4 = FiniteWall([w/2+h2, h/2+w2/2],[w/2+h2, h/2+r2], SV(1, 0))
ob3 = Semicircle([w/2+h2, h/2], r2, SV(-1, 0))
bd = Billiard(stembot, stemrightlow, stemrightupper,
 capbotright, cap, capbotleft, stemleft, ob1, ob5, ob2, ob4, ob3)

plot(bd)
# animate_evolution([randominside(bd) for i in 1:5], bd, 50)
# for i in 1:10
#     x, y = timeseries(randominside(bd), bd, 100)
#     plot(x, y)
# end
x, y = timeseries(randominside(bd), bd, 100)

# %% Make random particles:
ps = [Particle(0.0, h/2 + 0.01i, j) for i in -5:5 for j in 0:0.1:2π]
[plot(p) for p in ps]

κ = mean(parallelize(meancollisiontime, bd, 1000, ps))
