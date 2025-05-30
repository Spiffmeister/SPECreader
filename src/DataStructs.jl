

struct GeometryType{GT} end
const Toroidal = GeometryType{:Toroidal}()
const Cylindrical = GeometryType{:Cylindrical}()
const Slab = GeometryType{:Slab}()




struct VectorPotential{ATT}
    Ate :: ATT
    Aze :: ATT
    Ato :: ATT
    Azo :: ATT
end

struct FourierSeries{TT}
    R :: Vector{TT}
    Z :: Vector{TT}
end



"""
    SPECEquilibrium{TT, ATT<:Array{TT}, TSA}

All components needed to reconstruct a SPEC equilibrium.

```julia
eqname = "testing/data/G3V01L0Fi.002.sp.h5"
speceq = SPECEquilibrium(eqname)
```


"""
struct SPECEquilibrium{TT,ATT<:Array{TT},TSA}

    Ntor::Int64
    Geometry::GeometryType
    Mvol::Int64
    PoloidalResolution::Int64
    Nfp::Int64
    mn::Int64

    NumberofVolumes::Int64

    StellaratorSymmetric::Bool

    CoordinateSingularity::Bool

    m::Vector{Int64}
    n::Vector{Int64}

    Rbc::ATT
    Rbs::ATT
    Zbc::ATT
    Zbs::ATT

    Rpol::TT
    Rtor::TT

    RadialResolution::Vector{Int64}

    # Physics params

    Ate::ATT
    Aze::ATT
    Ato::ATT
    Azo::ATT


    # Storage for basis functions

    RadialBasis::TSA
    # Zernike2            :: Array{TT}
    # Zernike3            :: Array{TT}
    # Chebychev           :: Array{TT}

    cache::Vector{TT} # Temp storage to make methods non-allocating

    cache_metric::SArray{Tuple{3,3,3},TT}


    function SPECEquilibrium(eqname)

        f = h5open(eqname)

        foutput = f["output"]
        fphysics = f["input"]["physics"]



        # Geometry parameters

        Mvol = read(foutput["Mvol"])[1]

        mn = read(foutput["mn"])[1]

        m = read(foutput["im"])
        n = read(foutput["in"])

        # m .= convert.(Integer,m)


        Rbc = read(foutput["Rbc"])
        Zbs = read(foutput["Zbs"])
        Rbs = read(foutput["Rbs"])
        Zbc = read(foutput["Zbc"])


        TT = eltype(Rbc)

        # Physics parameters

        Igeometry = read(fphysics["Igeometry"])[1]

        if Igeometry == 3
            Geometry = Toroidal
        elseif Igeometry == 2
            Geometry = Cylindrical
        else
            Geometry = Slab
        end

        Geometry == Toroidal ? ICoordinateSingularity = true : ICoordinateSingularity = false

        StellaratorSymmetric = convert(Bool, read(fphysics["Istellsym"])[1])

        PoloidalResolution = read(fphysics["Mpol"])[1] # TODO : why is this here
        Ntor = read(fphysics["Ntor"])[1]
        Nfp = read(fphysics["Nfp"])[1]

        RadialResolution = read(fphysics["Lrad"])


        Nvol = read(fphysics["Nvol"])[1]

        Ate = read(f["vector_potential"]["Ate"])
        Aze = read(f["vector_potential"]["Aze"])
        Ato = read(f["vector_potential"]["Ato"])
        Azo = read(f["vector_potential"]["Azo"])



        Rpol = one(TT)
        Rtor = one(TT)

        try
            Rpol = read(fphysics["rpol"])[1]
            Rtor = read(fphysics["rtor"])[1]
        catch
        end

        RadialBasis = @MArray zeros(maximum(RadialResolution) + 1, PoloidalResolution + 1, 2)


        # n = Int.(n/Nfp)

        # VPList = (1,RadialResolution...)
        # VP = [VectorPotential(
        #         Ate[i:VPList[i+1],:],
        #         Aze[i:VPList[i+1],:],
        #         Ato[i:VPList[i+1],:],
        #         Azo[i:VPList[i+1],:]) for i in eachindex(VPList)]



        cache = zeros(TT, 3)

        new{TT,typeof(Ate),typeof(RadialBasis)}(Ntor, Geometry, Mvol, PoloidalResolution, Nfp, mn, Nvol,
            StellaratorSymmetric, ICoordinateSingularity,
            m, n,
            Rbc, Rbs, Zbc, Zbs, Rpol, Rtor,
            RadialResolution,
            Ate, Aze, Ato, Azo,
            RadialBasis,
            cache,
            @SArray zeros(TT, 3, 3, 3))

    end
end




"""
    ReadBoundary(fname::String)

Outputs the necessary data to reconstruct a SPEC boundary.
"""
function ReadBoundary(fname::String)

    f = h5open(fname)

    out = f["output"]

    phys = f["input"]["physics"]

    Rbc = read(out["Rbc"])
    Zbs = read(out["Zbs"])

    PoloidalResolution = read(phys["Mpol"])[end]
    Ntor = read(phys["Ntor"])[end]

    m = vcat([collect(-PoloidalResolution:PoloidalResolution) for i in 1:Ntor+1]...)[PoloidalResolution+1:end]
    n = vcat(zeros(Int,PoloidalResolution+1), repeat(1:Ntor,inner=2PoloidalResolution+1))

    if eltype(m) == Int32
        m = [promote(m...,Int64(1))[1:end-1]...]
    end

    boundary = (Rbc = Rbc[:,2],
        Zbs = Zbs[:,2],
        m = m,
        n = n)
    
    return boundary
end


"""
    ReadPoincare(fname::String)

Read the Poincaré section for a given SPEC equilibrium output.
"""
function ReadPoincare(fname::String)

    f = h5open(fname)

    pdata = f["poincare"]
    poinout = (R = read(pdata["R"]),Z = read(pdata["Z"]))

    return poinout
end
