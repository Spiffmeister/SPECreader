


"""
    choose_geometry
Specifies the geometry type based on input/physics/Igeometry
"""
function choose_geometry(geom)
    if geom == 1
        return Cartesian
    elseif geom == 2
        return Cylindrical
    elseif geom == 3
        return Toroidal
    else
        error("Unknown geometry")
    end
end

function readh5_single(Data,field)
    return read(Data,field)[1]
end


function read_spec(fname)

    fid = h5open(fname,"r")

    Inputs = fid["input/physics"]
    Outputs = fid["output"]

    # SPECequilibrium data struct
    Igeometry = choose_geometry(readh5_single(Inputs,"Igeometry"))
    # read(Inputs,"Istellsym")
    # read(Inputs,"Lfreebound")
    Nvol = readh5_single(Inputs,"Nvol")
    # read(Inputs,"Nfp")
    # read(Inputs,"Mpol")
    # read(Inputs,"Ntor")
    Lrad = read(Inputs,"Lrad") #Should always be a vector
    Mvol = readh5_single(Outputs,"Mvol")
    mn = readh5_single(Outputs,"mn")

    # Volume data struct
    Igeometry != Cartesian ? ICoordinateSingularity = false : ICoordinateSingularity = false

    Vol = Coordinates{Float64,Igeometry,Nvol}(
        read(Outputs,"Rbc"),
        read(Outputs,"Rbs"),
        read(Outputs,"Zbc"),
        read(Outputs,"Zbs"),
        ICoordinateSingularity,
        convert.(Float64,read(Outputs,"im")),
        convert.(Float64,read(Outputs,"in"))
    )


    Ate = read(fid["vector_potential"],"Ate")
    Aze = read(fid["vector_potential"],"Aze")
    Ato = read(fid["vector_potential"],"Ato")
    Azo = read(fid["vector_potential"],"Azo")

    SubVecs :: Array{VectorPotential} = []
    for ivol = 1:Mvol
        Ll = sum(Lrad[1:ivol-1]) + length(Lrad[1:ivol])
        Lu = sum(Lrad[1:ivol]) + length(Lrad[1:ivol])
        Vec = VectorPotential(
            Ate[Ll:Lu,1:mn],
            Aze[Ll:Lu,1:mn],
            Ato[Ll:Lu,1:mn],
            Azo[Ll:Lu,1:mn],
            Lrad[ivol],
            Igeometry,
            convert(typeof(Lrad[ivol]),ivol)
        )
        push!(SubVecs,Vec)
    end
    
    SE = SPECequilibrium(
        Igeometry, 
            convert(Bool,readh5_single(Inputs,"Istellsym")), convert(Bool,readh5_single(Inputs,"Lfreebound")), 
            readh5_single(Inputs,"Nvol"), 
            readh5_single(Inputs,"Nfp"), readh5_single(Inputs,"Mpol"), readh5_single(Inputs,"Ntor"), 
            Lrad, readh5_single(Outputs,"Mvol"), mn,
        SubVecs,
        Vol
    )

    close(fid)

    return SE
end