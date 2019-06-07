import FileIO, BSON, DrWatson, Pkg

Pkg.activate(dirname(@__DIR__))
@assert projectname() == "MagneticBilliardsLyapunovs"

function symboldict(d::Dict{String, <:Any})
    c = Dict{Symbol, Any}()
    for e in d
        c[Symbol(e[1])] = e[2]
    end
    c
end

function change(path)
    readdir(path)
    for f in readdir(path)
        if f[end-3:end] == "jld2"
            data = FileIO.load(joinpath(path, f))
            c = symboldict(data)
            bsonname = joinpath(path, f[1:end-4]*"bson")
            BSON.bson(bsonname, c)
            # rm(joinpath(path, f))
        end
    end
end

change(datadir()*"mushrooms")
