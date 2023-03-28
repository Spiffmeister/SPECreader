
"""
    Geometry
Cartesian, Cylindrical or Toroidal
"""
struct GeometryType{G} end
const Cartesian = GeometryType{:CartesianGeometry}()
const Cylindrical = GeometryType{:CylindricalGeometry}()
const Toroidal = GeometryType{:ToroidalGeometry}()

"""
    Parameters
Stores simulation output parameters
"""
struct Parameters
    # Input-list physics
    StellaratorSymmetry    :: Bool
    Lfreebound              :: Bool
    Nvol                    :: Integer
    Nfp                     :: Integer
    Mpol                    :: Integer
    Ntor                    :: Integer
    Lrad                    :: AbstractArray{Integer}
    # Output-list
    Mvol                    :: Integer
    mn                      :: Integer
end


"""
    VectorPotential
Stores the data needed to construct the vector potential in a given volume
"""
struct VectorPotential{T,N,NV}

    Ate         :: AbstractArray{T}
    Ato         :: AbstractArray{T}
    Aze         :: AbstractArray{T}
    Azo         :: AbstractArray{T}

    Lrad        :: N

    Isingular   :: Bool # Singularity present (inner most volume)

    function VectorPotential(Ate::AbstractArray{T},Ato::AbstractArray{T},Aze::AbstractArray{T},Azo::AbstractArray{T},
            Lrad::N,Geom::GeometryType,VolId::N) where {T,N}
        new{T,N,VolId}(
            Ate,Ato,Aze,Azo,
            Lrad,
            (VolId == 1) & (Geom != Cartesian) ? true : false
        )
    end
end

"""
    Metrics
Storage for metric computation
"""
mutable struct Metrics{T}
    x       :: AbstractArray{T}

    Jac     :: T
    ∇Jac    :: AbstractArray{T}
    JacMat  :: AbstractArray{T}

    gᵢⱼ     :: AbstractArray{T}
    dgᵢⱼ    :: AbstractArray{T}

    R₀₀     :: T
    Rᵢⱼ     :: AbstractArray{T}
    Zᵢⱼ     :: AbstractArray{T}

    function Metrics(T)
        x       = zeros(T,3)
        Jac     = T(0) 
        ∇Jac    = zeros(T,(3,3))
        JacMat  = zeros(T,(3,3))
        gᵢⱼ     = zeros(T,(3,3))
        dgᵢⱼ    = zeros(T,(3,3,3))
        R₀₀     = T(0)
        Rᵢⱼ     = zeros(T,(4,3))
        Zᵢⱼ     = zeros(T,(4,3))

        new{T}(x,
            Jac, ∇Jac, JacMat,
            gᵢⱼ, dgᵢⱼ,
            R₀₀, Rᵢⱼ, Zᵢⱼ)
    end
end

"""
    Coordinates
Data for constructing coordinates
"""
struct Coordinates{T,G,V}
    Rbc     :: AbstractArray{T}
    Rbs     :: AbstractArray{T}
    Zbc     :: AbstractArray{T}
    Zbs     :: AbstractArray{T}

    CoordinateSingularity   :: Bool
    mᵢ                      :: AbstractArray{T}
    nᵢ                      :: AbstractArray{T}

end

"""
    SPECequilibrium
SPEC equilibrium hdf5 data
"""
struct SPECequilibrium{G,V}
    Params  :: Parameters
    A       :: Vector{VectorPotential}
    Ri      :: Coordinates
    M       :: Metrics

    function SPECequilibrium(Igeometry::GeometryType,StellSym::Bool,Freebound::Bool,
        Nvol::N,Nfp::N,Mpol::N,Ntor::N,Lrad::AbstractArray{N},
        Mvol::N,mn::N,
        A::Vector{VectorPotential},
        Vol::Coordinates) where {N}

        Params = Parameters(StellSym::Bool,Freebound::Bool,
        Nvol::N,Nfp::N,Mpol::N,Ntor::N,Lrad::AbstractArray{N},
        Mvol::N,mn::N)

        Met = Metrics(typeof(A[1].Ate[1]))

        new{Igeometry,Nvol}(
            Params,
            A,
            Vol,
            Met)
    end
end

"""
    GetType
"""
GetType(SE::SPECequilibrium{G,V}) where {G,V} = G

"""
    GetVol
"""
function GetVol end
GetVol(A::VectorPotential{T,N,V}) where {T,N,V} = V
GetVol(SE::SPECequilibrium{G,V}) where {G,V} = V




Base.show(io::IO, SE::SPECequilibrium) = print(io,GetType(SE)," SPEC equilibrium read with ",GetVol(SE)," volumes.")