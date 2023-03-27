
"""
    GetMetric
Get the metric information at (s,θ,ξ)
"""
function GetMetric end
function GetMetric(M::Metrics,s::T,θ::T,ξ::T,vᵢ::Integer,Coords::Coordinates{T,Cartesian,V},P::Parameters) where {T,V}

    R_Derivative!(M,s,θ,ξ,vᵢ,Coords,P)

    M.gᵢⱼ[:,:]  = T(0)
    M.J         = T(0)
    M.dgᵢⱼ[:,:] = T(0)

    M.Jac   = Rᵢⱼ[1,1]
    M.∇Jac  = Rᵢⱼ[1:3,1]

    M.JacMat[1,2]   = T(1)
    M.JacMat[2,3]   = T(1)
    M.JacMat[3,1:3] = Rᵢⱼ[1,1:3]

    for i = 1:3
        for j = 1:3
            gᵢⱼ[j,i] += Rᵢⱼ[1,j]*Rᵢⱼ[1,i]
            for k = 1:3
                dgᵢⱼ[i,j,k] += Rᵢⱼ[k,j]*Rᵢⱼ[1,i] + Rᵢⱼ[1,j]*Rᵢⱼ[k,i]
            end
        end
    end

    gᵢⱼ[2,2] += T(1)
    gᵢⱼ[3,3] += T(1)

    return 
end

"""
    R_Derivative
Compute R and its derivatives
"""
function R_Derivative!(M::Metrics,s::T,θ::T,ξ::T,vᵢ::Integer,Coord::Coordinates,P::Parameters) where {T}

    a      = αᵢ(Coord.mᵢ,Coord.nᵢ,θ,ξ)
    cosα   = cos.(a)
    sinα   = sin.(a)

    t = zeros(T,(4,P.mn))

    alss = (1.0-s)/2.0
    blss = (1.0+s)/2.0

    t[1,:] = alss*Coord.Rbc[:,vᵢ] + blss*Coord.Rbc[:,vᵢ+1]
    t[2,:] = -0.5*Coord.Rbc[:,vᵢ] + 0.5*Coord.Rbc[:,vᵢ+1]

    if !P.Istellaratorsymmetry
        t[3,:] = alss*Coord.Rbs[:,vᵢ] + blss*Coord.Rbs[:,vᵢ+1]
        t[4,:] = -0.5*Coord.Rbs[:,vᵢ] + 0.5*Coord.Rbs[:,vᵢ+1]
    end


    M.R₀₀       = dot(t[1,:],cosa)                      #R coordinate
    
    M.Rᵢⱼ[1,1]  = dot(t[2,:],cosa)                      #dRds
    M.Rᵢⱼ[1,2]  = dot(t[1,:],-Coord.mᵢ*sina)            #dRdu
    M.Rᵢⱼ[1,3]  = dot(t[1,:],Coord.nᵢ*sina)             #dRdv
    
    M.Rᵢⱼ[2,1]  = T(0)                                  #d2Rds2
    M.Rᵢⱼ[2,2]  = dot(t[2,:],-Coord.mᵢ.*sina)           #d2Rdsdu
    M.Rᵢⱼ[2,3]  = dot(t[2,:],Coord.nᵢ.*sina)            #d2Rdsdv
    
    M.Rᵢⱼ[3,2]  = dot(t[1,:],-Coord.mᵢ.^2 .*cosa)       #d2Rdu2
    M.Rᵢⱼ[3,3]  = dot(t[1,:],Coord.mᵢ.*Coord.nᵢ .*cosa) #d2Rdudv

    M.Rᵢⱼ[4,3]  = dot(t[1,:],-Coord.nᵢ.^2 .*cosa)       #d2Rdv2

    if !P.Stellaratorsymmetry
        M.R₀₀       += dot(t[3,:],sina)

        M.Rᵢⱼ[1,1]  += dot(t[4,:],sina)
        M.Rᵢⱼ[1,2]  += dot(t[3,:],Coord.mᵢ.*cosa)
        M.Rᵢⱼ[1,3]  += dot(t[3,:],-Coord.nᵢ.*cosa)

        M.Rᵢⱼ[2,1]  += dot(t[3,:],sina)
        M.Rᵢⱼ[2,2]  += dot(t[3,:],sina)
        M.Rᵢⱼ[2,3]  += dot(t[3,:],sina)

        M.Rᵢⱼ[3,2]  += dot(t[3,:],sina)
        M.Rᵢⱼ[3,3]  += dot(t[3,:],sina)
        M.Rᵢⱼ[4,3]  += dot(t[3,:],sina)
    end


    # Rᵢⱼᵀ = Rᵢⱼ
    M.Rᵢⱼ[3,1] = M.Rᵢⱼ[2,2]
    M.Rᵢⱼ[4,1] = M.Rᵢⱼ[2,3]
    M.Rᵢⱼ[4,2] = M.Rᵢⱼ[3,3]

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