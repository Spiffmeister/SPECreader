
"""
    Geometry
"""
struct GeometryType{G} end
const Cartesian = GeometryType{:CartesianGeometry}()
const Cylindrical = GeometryType{:CylindricalGeometry}()
const Toroidal = GeometryType{:ToroidalGeometry}()



"""
    VectorPotential
"""
struct VectorPotential{T,N}

    Ate         :: AbstractArray{T}
    Ato         :: AbstractArray{T}
    Aze         :: AbstractArray{T}
    Azo         :: AbstractArray{T}

    Lrad        :: N

    isingular   :: Bool # Singularity present (inner most volume)
end

"""
    Metric

Holds Jacobian
"""
struct Metric{T}
    Jac     :: T
    x       :: AbstractArray{T}
    ∇Jac    :: AbstractArray{T}
    gᵢⱼ     :: AbstractArray{T}
    dgᵢⱼ    :: AbstractArray{T}

    Rᵢⱼ     :: AbstractArray{T}
    Zᵢⱼ     :: AbstractArray{T}
end

"""
    Volume
"""
struct Geom{T,G,V}
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
T <: Type, V <: Number of volumes
"""
struct SPECequilibrium{G,V}
    # Input-list physics
    Igeometry               :: GeometryType
    Istellaratorsymmetry    :: Bool
    Lfreebound              :: Bool
    Nvol                    :: Integer
    Nfp                     :: Integer
    Mpol                    :: Integer
    Ntor                    :: Integer
    Lrad                    :: AbstractArray{Integer}
    # Output-list
    Mvol                    :: Integer
    mn                      :: Integer

    A                       :: Vector{VectorPotential}
    Ri                      :: Geom

    function SPECequilibrium(Igeometry::GeometryType,StellSym::Bool,Freebound::Bool,
        Nvol::N,Nfp::N,Mpol::N,Ntor::N,Lrad::AbstractArray{N},
        Mvol::N,mn::N,
        A::Vector{VectorPotential},
        Vol::Geom) where {N}


        new{Igeometry,Nvol}(Igeometry,
            StellSym,Freebound,
            Nvol,Nfp,Mpol,Ntor,Lrad,
            Mvol,mn,
            A,
            Vol)
    end
end


GetType(SE::SPECequilibrium{G,V}) where {G,V} = G



Base.show(io::IO, SE::SPECequilibrium) = print(io,GetType(SE)," SPEC equilibrium read with ",SE.Nvol," volumes.")