


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

function readh5(Data,field)
    if size(Data[field]) == (1,)
        return read(Data,field)[1]
    else
        return read(Data,field)
    end
end


function read_spec(fname)

    fid = h5open(fname,"r")

    Inputs = fid["input/physics"]
    Outputs = fid["output"]

    # SPECequilibrium data struct
    Igeometry = choose_geometry(readh5(Inputs,"Igeometry")[1])
    # read(Inputs,"Istellsym")
    # read(Inputs,"Lfreebound")
    Nvol = readh5(Inputs,"Nvol")
    # read(Inputs,"Nfp")
    # read(Inputs,"Mpol")
    # read(Inputs,"Ntor")
    Lrad = read(Inputs,"Lrad") #Should always be a vector
    Mvol = readh5(Outputs,"Mvol")
    mn = readh5(Outputs,"mn")

    # Volume data struct
    Igeometry != Cartesian ? ICoordinateSingularity = false : ICoordinateSingularity = false
    
    Vol = Geom{Float64,Igeometry,Nvol}(
        readh5(Outputs,"Rbc"),
        readh5(Outputs,"Rbs"),
        readh5(Outputs,"Zbc"),
        readh5(Outputs,"Zbs"),
        ICoordinateSingularity,
        readh5(Outputs,"im"),
        readh5(Outputs,"in")
    )


    Ate = readh5(fid["vector_potential"],"Ate")
    Aze = readh5(fid["vector_potential"],"Aze")
    Ato = readh5(fid["vector_potential"],"Ato")
    Azo = readh5(fid["vector_potential"],"Azo")

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
            (ivol == 1) & (Igeometry != Cartesian) ? true : false)
        push!(SubVecs,Vec)
    end
    
    SE = SPECequilibrium(
        Igeometry, 
            convert(Bool,readh5(Inputs,"Istellsym")), convert(Bool,readh5(Inputs,"Lfreebound")), 
            readh5(Inputs,"Nvol"), 
            readh5(Inputs,"Nfp"), readh5(Inputs,"Mpol"), readh5(Inputs,"Ntor"), 
            Lrad, readh5(Outputs,"Mvol"), mn,
        SubVecs,
        Vol
    )

    close(fid)

    return SE
end