using BSON, PyPlot

# this is my formula for the phasespace volume
ps_volume(l = 1.0, w = 1.0, ω = 1.0; B  = ω/2) = 2π * ((1+π)/B^2 + ((1+2π)*w + 2l)/B + w*l + π*w^2)

"""
    affine(d1, d2)
Scales d1 to match d2. Doesn't shift the data so we can actually see the shift.
"""
function affine(d1, d2)
    a,b = extrema(d1)
    c,d = extrema(d2)

    n = (a*d - c*b)/(a - b)
    m = (c - n)/a

    println("Transform: $(m)*x + $(n)")
    return m.*d1 #.+ n
end

function make_plot(data; complicated = false)
    # make sure that everything is sorted by ω
    sp = sortperm(data[:ωs])

    #convenience
    ωs = data[:ωs][sp]
    λss = data[:λss][sp]
    σλs = data[:σλs][sp]
       
    # simple method of getting volume - just don't normalise with total
    vc_simple = (data[:psv_vis][sp])

    # use num. g_c and multiply analyt. V_tot
    pv = ps_volume.(1.0, 1.0, ωs)
    vc_complicated = (data[:psv_vis][sp]./data[:psv_tot][sp]) .* pv

    f1 = figure("lyaps")
    ax = plt.gca() #plt.twinx()

    # plot Lyapunovs with errorbar
    ax.errorbar(ωs, λss, σλs, marker = "s", color = "C0", label=L"\lambda")

    # plot one over V_c
    ax.plot(ωs[sp], affine(1 ./vc_simple, λss), "o", color = "C1", label=L"$1/V_C$ from unnormalised box-counting", zorder = 10)
    complicated && ax.plot(ωs[sp], affine(1 ./vc_complicated, λss), "^", color = "C2", label=L"$1/V_C$ from normalised box-counting & analyt. $V_{tot}$", zorder = 20)

    
    ax.set_xlabel(L"\omega")
    ax.set_ylabel(L"1/t")
    ax.semilogx()

    ax.fill_between(ωs, λss[1]*ones(size(λss)), λss, color = "C0", alpha = 0.1)
    plt.annotate(s="", xy=(ωs[end],λss[1] ), xytext=(ωs[end],λss[end]), arrowprops=Dict(:arrowstyle=>"<->", :lw=>2, :shrinkB=>0))
    plt.text(s=L"\delta\lambda", ωs[end]*0.6  , (λss[end]-λss[1])/2, fontsize = 15)
    f1.tight_layout()

    f1.legend(loc = (0.3, 0.2))
end

data_normal = BSON.load("output_Tbox=100.0_Nbox=100.bson")
data_better = BSON.load("output_Tbox=10000.0_Nbox=200.bson")
data_best = BSON.load("output_Tbox=100000.0_Nbox=300.bson")


make_plot(data_best; complicated=true)
#=
f2 = figure("volumes")
plt.plot(ωs[sp], psts[sp]./psts[end] ,"o",  label="from num. normalisation factor", color = "C0")
plt.plot(ωs[sp], pv./pv[end], label=L"from analytic $V_{tot}$", color = "C1")
axhline(1.0, label=L"Limit for $ω \rightarrow \infty$", color = "C2")
xlabel(L"\omega")
ylabel(L"V_{tot}(ω)/V_{tot}(10^3)")
title("Total phase space volume")
semilogx()
#f1.legend(loc = (0.45, 0.125))
tight_layout()

ins1 = f2.add_axes([0.5, 0.5, 0.4, 0.36])
ins1.set_title("rescaled y axes")
ins1.plot(ωs[sp], psts[sp]./psts[end] ,"o",  label="from normalisation factor", color = "C0")
ylabel("num.")
xlabel(L"\omega")
semilogx()
ins2 = ins1.twinx()
ylabel("analyt.")
ins2.plot(ωs[sp], pv./pv[end], label=L"from analytic $V_c$", color = "C1")
=#
