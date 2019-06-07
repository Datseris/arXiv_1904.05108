using DrWatson
quickactivate(@__DIR__, "MagneticBilliardsLyapunovs")
using Distributed

nworkers() == 1 && addprocs()
@everywhere using DynamicalBilliards

@everywhere function mean_cell_times(h = 1.0, w = 0.1; T = 10000.0, N = 1000)

    bd = billiard_mushroom(h, w, 1.0)
    capobsts = (3, 4, 5); cap = 4
    stemobsts = (1, 2, 6)

    tcs, tls, tlcs = 0.0, 0.0, 0.0
    nl, nc = 0, 0 # number of cells

    for j ∈ 1:N # particle loop
        p = MushroomTools.randomchaotic(h, w, 1.0)

        # bring particle to cap before counting time
        obst, = bounce!(p, bd)
        while obst != cap
            obst, = bounce!(p, bd)
        end
        tt = zero(T)

        while tt < T
            t_ret = 0.0
            t_help = 0.0
            # Bounce once to see what cell we are in:
            # (previous obstacle was definitely a cap head)
            obst, t, = bounce!(p, bd)
            tt += DynamicalBilliards.increment_counter(tt, t)
            t_ret = t
            if obst ∈ capobsts # we are in a laminar cell
                while obst ∉ stemobsts
                    obst, t, = bounce!(p, bd)
                    tt += DynamicalBilliards.increment_counter(tt, t)
                    t_ret += t
                    t_help = t
                end
                # notice laminar cell does not include last collision
                # (because it is not in cap)
                tls += (t_ret - t_help)
                t_ret = t_help
                while obst ≠ cap
                    obst, t, = bounce!(p, bd)
                    tt += DynamicalBilliards.increment_counter(tt, t)
                    t_ret += t
                end
                tlcs += t_ret
                nl += 1
            else # we are in a chaotic cell:
                # bounce back to cap and add times
                while obst ≠ cap
                    obst, t, = bounce!(p, bd)
                    tt += DynamicalBilliards.increment_counter(tt, t)
                    t_ret += t
                end
                tcs += t_ret
                nc += 1
            end
        end
    end

    println("finished w = $(w)")
    # return total/# data points
    return (tcs/nc, tls/nl, tlcs/nl, nl/(nl+nc))
end

T = 20000.0; N = 20000

ws = 0.1:0.1:2.0
const hhh = 1.0
ret = @time pmap(w->mean_cell_times(hhh, w; T = T, N = N), ws)

a = [zeros(length(ws)) for i in 1:length(ret[1])]
for i in 1:length(ret[1])
    a[i] = [x[i] for x in ret]
end
tc, tl, tlc, fl = a


# %% plot
# Geometric return time to stem, without laminar phases:
ταυ_c(h, w) =  2*(1 - w^2/(36) - w^4/(1200) - w^6/(15680) - w^8/(145152)) + π*h
# mean return time to stem, from Kacs lemma, includes laminar phases
καππα(h, w) =  (π/2*(2h + π/2)) * 2MushroomTools.V_3D_cha(h, w, 1.0) / (2π*w*(2h + π/2))

using PyPlot
figure(figsize = (16, 8))
subplot(1,2,1)
scatter(ws, tc, color = "C0", label = "<tc>")
scatter(ws, @.((1-fl)*tc + fl*tlc), color = "C0", marker ="D", label = "(1-p)<tc> + (p)<tlc>")
plot(ws, ταυ_c.(hhh, ws), ls = "dashed", label = "\$t_c\$", color = "C0")
scatter(ws, @.((1-fl)*tc + fl*(tl + tlc)), color = "C1", marker ="D", label = "(1-p)<tc> + (p)<tl + tlc>")
legend()


subplot(1,2,2)
plot(ws, καππα.(hhh, ws), ls = "dashed", label = "\$\\kappa\$", color = "C1")
# scatter(ws, @.((1 - fl)*tc + fl*tll), label = "(1-p)<tc> + p(<tl> + <tlc>)",
# marker = "D", color = "C1")
scatter(ws, @.((1-fl)*tc + fl*tl), label = "(1-p)<tc> + p<tl>",
marker = "s", color = "C2")
scatter(ws, @.(tc + fl*tl), label = "<tc> + p<tl>",
marker = "X", color = "C3")
scatter(ws, @.(tc + fl*(tl + tlc)), label = "<tc> + p(<tl> + <tlc>)",
marker = "h", color = "C4")
scatter(ws, @.((1-fl)*tc + fl*(tl + tlc)),
label = "(1-p)<tc> + p(<tl> + <tlc>)",
marker = "o", color = "C5")

legend()
