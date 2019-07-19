################################################################################
##
##  This file should produce a plot similar to Fig. 3.8 from the thesis.
##
################################################################################

using DrWatson
quickactivate(@__DIR__, "MagneticBilliardsLyapunovs")
using Distributed

# add worker processes if they don't already exist
nworkers() == 1 && addprocs()

# definitions needed on all workers
@everywhere using DynamicalBilliards

@everywhere function get_stemreturn(l = 1.0, w = 0.1, r = 1.0; T = 10000.0, N = 10000)
    bd = billiard_mushroom(l, w, r)

    # save return times here
    rettimes = 0.0
    # count no. of data points
    rets = 0

    for j ∈ 1:N
        # create chaotic particle
        p = MushroomTools.randomchaotic(l, w, r)

        # bring particle to stem bottom
        obst, = bounce!(p, bd)
        while obst != 1
            obst, = bounce!(p, bd)
        end

        tt = 0.0
        while tt < T
            t_ret = 0.0
            obst = 0
            # bounce until back at stem bottom & keep track of time
            while obst != 1
                obst, t, = bounce!(p, bd)
                tt += t
                t_ret += t
            end
            # add stem return time to total
            rettimes += t_ret
            rets += 1
        end
    end
    # status report
    println("finished w = $(w)")
    # return total/# data points
    return rettimes/(rets)
end

T = 10000.0; N = 10000

ws = 0.1:0.05:2.0
κs = @time pmap(w->get_stemreturn(1.0, w, 1.0; T=T, N=N), ws)

# %% plot
V_3D_ch(l,w,r) = MushroomTools.V_3D_tot(l,w,r) - MushroomTools.V_3D_reg(l,w,r)
# note the factor of two missing from my thesis
κs_kac(l,w,r) = (π/2*(2l + π/2)) * 2*V_3D_ch(l,w,r)/(2π * w*r*(2l + π*r/2))

figure()
xlabel(L"w")
ylabel(L"κ_s")

plot(ws, κs, "o", label="numerical data")
plot(ws, κs_kac.(1.0, ws, 1.0), label="Kac's lemma")
text(0.1, 0.5, "h=1.0, T=$T, N=$N", transform = gca().transAxes)
title("thesis fig 3.8: mean return time to stem")
tight_layout()

legend()
savefig(plotsdir()*"mushrooms/fig38.png")

# remove additional processes
# rmprocs(workers())
