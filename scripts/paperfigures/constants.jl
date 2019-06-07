using PyPlot
SUBPLOT = 24
mixedls = [ ":", "-", "--", "-."]
markers = ["o", "s", "^", "p", "P", "D"]
PyPlot.rc("lines", lw = 2.5)
PyPlot.rc("font", size = 24) # set default fontsize
PyPlot.rc("legend", fontsize = 24)
const figx = 12 # width of all figures

mutable struct CoolColors
    c::Vector
    n::Int
end
CoolColors(c) = CoolColors(c, 0)

Base.getindex(c::CoolColors, i) = c.c[(i-1)%length(c.c) + 1]
function Base.getindex(c::CoolColors)
    c.n += 1
    c[c.n]
end
coolcolors = CoolColors(
["#65ADC2", "#233B43", "#E84646",
"#C29365", "#168E7F", "#985CC9",
 "#822e01", "#7c6d8c"],
0)


# figure()
# for i in 1:length(coolcolors.c)
#     plot([0, 1], [0, 1] .+ i, color = coolcolors[i], label = "$i")
# end
# legend()

function add_identifiers!(fig = gcf())
    bbox = Dict(:boxstyle => "round,pad=0.3", :facecolor=>"white", :alpha => 1.0)
    for (i, ax) in enumerate(fig.get_axes())
        l = collect('a':'z')[i]
        x = 0.985
        loc = (x, 0.975)
        ax.text(loc..., "$(l)", size = SUBPLOT,
        transform = ax.transAxes, bbox = bbox, zorder = 99)
    end
end

function coolfill(x, y, dy, ax, c, label = "")
    α = 0.25
    ax.plot(x, y, label = label, color = c, lw = 2.0)
    ax.fill_between(x, y .- dy, y .+ dy, alpha = α, color = c)
    lw = 0.5
    α2 = 0.5
    ax.plot(x, y .+ dy,  color = c, lw = lw, alpha = α2)
    ax.plot(x, y .- dy,  color = c, lw = lw, alpha = α2)
end

function nice_arrow(xc, yc, xspan, yspan, ax = gca();
    style = "<->", tex = "", xo = 0.2, yo = -0.2)
    ax.annotate("",  xy=(xc-xspan/2, yc - yspan/2), xycoords="data",
                xytext=(xc+xspan/2, yc + yspan/2), textcoords="data",
                arrowprops = (Dict(:arrowstyle=>style,
                                :connectionstyle=>"arc3",
                                :lw=>1.5, :facecolor => "black")), zorder = 99)
    if tex != ""
        ax.text(xc + xo, yc + yo, tex, size = 24)
    end
end

function coolhist(ax, data, bins, color, label = "", alpha = 0.25)
    h, b, = ax.hist(data, bins, density = true, color = color,
    alpha = alpha)

    b = 0.5(b[1:end-1] .+ b[2:end])
    ax.plot(b, h, color = color, lw = 1.0, label = label)
end

function add_grid!(ax, n::Int; kwargs...)
    @assert n ≥ 3
    x = ax.get_xlim()
    y = ax.get_ylim()
    dx = (x[2]-x[1])/n; dy = (y[2]-y[1])/n
    for i in 1:(n-1)
        ax.axhline(y[1]+n*dy, color = "gray", alpha = 0.5)
        ax.axvline(x[1]+n*dx, color = "gray", alpha = 0.5)
    end
end
