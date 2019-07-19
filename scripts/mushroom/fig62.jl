################################################################################
##
##  This file should produce a plot similar to Fig. 6.2 from the thesis.
##
################################################################################

using DrWatson
quickactivate(@__DIR__, "MagneticBilliardsLyapunovs")
using Distributed

# add worker processes if they don't already exist
nworkers() == 1 && addprocs()

# definitions needed on all workers
@everywhere using DynamicalBilliards

# Importantly this function includes both types of chaotic cells
@everywhere function get_capviastem(l = 1.0, w = 0.1, r = 1.0; T = 10000.0, N = 1000)
    bd = billiard_mushroom(l, w, r)

    # save return times here
    rettimes = 0.0
    # count no. of data points
    rets = 0

    for j ∈ 1:N
        # create chaotic particle
        p = MushroomTools.randomchaotic(l, w, r)

        # bring particle to cap
        obst, = bounce!(p, bd)
        while obst != 4
            obst, = bounce!(p, bd)
        end

        tt = 0.0
        while tt < T
            t_ret = 0.0
            # bounce until leaving cap & keep track of last collision time
            while obst ∈ (3,4,5)
                obst, t, = bounce!(p, bd)
                tt += t
                t_ret = t
            end

            # bounce back to cap and add times
            while obst != 4
                obst, t, = bounce!(p, bd)
                tt += t
                t_ret += t
            end

            # add return time to total
            rettimes += t_ret
            rets += 1
        end
    end
    # status report
    println("finished w = $(w)")
    # return total/# data points
    return rettimes/rets
end

T = 10000.0; N = 10000

ws = 0.1:0.05:2.0
tcs = @time pmap(w->get_capviastem(1.0, w, 1.0; T = T, N = N), ws)

# %%
using PyPlot
figure()
xlabel(L"w")
ylabel(L"t_c")

plot(ws, tcs, "o", label="numerical data")
plot(ws, fill(2 + π, length(ws)), label="very simple approximation")
plot(ws, 2 .* (@. 1 - ws^2/(36) - ws^4/(1200) - ws^6/(15680) - ws^8/(145152)) .+ π,
     label="simple approximation")

legend()
text(0.1, 0.5, "h=1.0, T=$T, N=$N", transform = gca().transAxes)
title("thesis fig 6.2: mean chaotic cell time")
tight_layout()
savefig(plotsdir()*"mushrooms/fig62.png")
# remove additional processes
# rmprocs(workers())
