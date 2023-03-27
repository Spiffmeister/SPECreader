
"""
    GetMetric
Get the metric information at (s,θ,ξ)
"""
function GetMetric end
function GetMetric(s::T,θ::T,ξ::T,Coords::Coordinates{T,G::GeometryType{:CartesianGeometry},V},P::Parameters) where {T,V}

    R₀₀,Rᵢⱼ = R_Derivative(s,θ,ξ,vᵢ,Coords,P)

    gᵢⱼ     = T(0)
    J       = T(0)
    dgᵢⱼ    = T(0)

    Jac = Metric.Rᵢⱼ[1,2]
    ∇Jac = Metric.Rᵢⱼ[2:4,2]
    
    J = zeros(T,(3,3))
    J[1,1] = J[2,2] = 1
    J[3,1:3] = Metric.Rᵢⱼ[1,2:4]

    for i = 1:3
        for j = 1:3
            gᵢⱼ[j,i] += Metric.Rᵢⱼ[1,j]*Metric.Rᵢⱼ[1,i]
            for k = 1:3
                dgᵢⱼ[i,j,k] += Metric.Rᵢⱼ[k,j]*Metric.Rᵢⱼ[1,i] + Metric.Rᵢⱼ[1,j]*Metric.Rᵢⱼ[k,i]
            end
        end
    end

    gᵢⱼ[3,3] += 1.
    gᵢⱼ[4,4] += 1.

end

"""
    R_Derivative
Compute R and its derivatives
"""
function R_Derivative(s::T,θ::T,ξ::T,vᵢ::Integer,Coord::Coordinates{T,G,V},P::Parameters) where {T}

    a      = αᵢ(Coord.mᵢ,Coord.nᵢ,θ,ξ)
    cosα   = cos.(a)
    sinα   = sin.(a)

    t = zeros(T,(4,P.mn))
    ddt = zeros(T,(2,P.mn))


    if Coord.CoordinateSingularity & (lvol==1)
        sbar = (1.0 + s) / 2.0


    else

        alss = (1.0-s)/2.0
        blss = (1.0+s)/2.0

        t[1,:] = alss*Coord.Rbc[:,vᵢ-1] + blss*Coord.Rbc[:,vᵢ]
        t[2,:] = -0.5*Coord.Rbc[:,vᵢ-1] + 0.5*Coord.Rbc[:,vᵢ]

        # ddt1 = 0.0

        if !P.Istellaratorsymmetry
            t[3,:] = alss*Coord.Rbs[:,vᵢ-1] + blss*Coord.Rbs[:,vᵢ]
            t[4,:] = -0.5*Coord.Rbs[:,vᵢ-1] + 0.5*Coord.Rbs[:,vᵢ]
        end

    end

    R₀₀ = dot(t[1,:],cosa) #R coordinate
    
    Rᵢⱼ[1,1] = dot(t[2,:],cosa)             #dRds
    Rᵢⱼ[1,2] = dot(t[1,:],-Coord.mᵢ*sina)   #dRdu
    Rᵢⱼ[1,3] = dot(t[1,:],Coord.nᵢ*sina)    #dRdv
    
    Rᵢⱼ[2,1] = T(0)                         #d2Rds2
    Rᵢⱼ[2,2] = dot(t[2,:],-Coord.mᵢ.*sina)  #d2Rdsdu
    Rᵢⱼ[2,3] = dot(t[2,:],Coord.nᵢ.*sina)   #d2Rdsdv
    
    Rᵢⱼ[3,1] = dot(t[1,:],-Coord.mᵢ.^2 .*cosa)          #d2Rdu2
    Rᵢⱼ[3,2] = dot(t[1,:],Coord.mᵢ.*Coord.nᵢ .*cosa)    #d2Rdudv
    Rᵢⱼ[3,3] = dot(t[1,:],-Coord.nᵢ.^2 .*cosa)          #d2Rdv2

    if !P.Stellaratorsymmetry
        R₀₀ += dot(t[3,:],sina)

        Rᵢⱼ[1,1] += dot(t[4,:],sina)
        Rᵢⱼ[1,2] += dot(t[3,:],Coord.mᵢ.*cosa)
        Rᵢⱼ[1,3] += dot(t[3,:],-Coord.nᵢ.*cosa)

        Rᵢⱼ[2,1] += dot(t[3,:],sina)
        Rᵢⱼ[2,2] += dot(t[3,:],sina)
        Rᵢⱼ[2,3] += dot(t[3,:],sina)

        Rᵢⱼ[2] += dot(t[3,:],sina)
        Rᵢⱼ[2] += dot(t[3,:],sina)
        Rᵢⱼ[2] += dot(t[3,:],sina)
    end


    # Rᵢⱼᵀ = Rᵢⱼ
    Rᵢⱼ[2,1] = Rᵢⱼ[1,2]
    Rᵢⱼ[3,1] = Rᵢⱼ[1,3]
    Rᵢⱼ[3,2] = Rᵢⱼ[2,3]

    
    return R₀₀,Rᵢⱼ
end

"""
    Z_Derivative
"""
function Z_Derivative()

    # cosa = cos.(αᵢ())
    # sina = sin.(αᵢ())

end


@inline αᵢ(m,n,θ,ξ) = m*θ - n*ξ




function Rᵢⱼ(R,s,θ,ζ,
    mᵢ,nᵢ,Rbs,Rbc,Ntor,vᵢ)

end