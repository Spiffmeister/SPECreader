
"""
    GetMetric
Get the metric information at (s,θ,ξ)
"""
function GetMetric end
function GetMetric(SE::SPECequilibrium{Cartesian,V},s::T,θ::T,ξ::T,vᵢ::Integer) where {T,V}

    R_Derivative!(SE.M,s,θ,ξ,vᵢ,SE.Ri,SE.Params)

    SE.M.x[1] = θ
    SE.M.x[2] = ξ
    SE.M.x[3] = SE.M.R₀₀

    SE.M.gᵢⱼ[:,:]  .= T(0)
    SE.M.dgᵢⱼ[:,:,:] .= T(0)

    SE.M.Jac   = SE.M.Rᵢⱼ[1,1]
    SE.M.∇Jac  = SE.M.Rᵢⱼ[1:3,1]

    SE.M.JacMat[1,2]   = T(1)
    SE.M.JacMat[2,3]   = T(1)
    SE.M.JacMat[3,1:3] = SE.M.Rᵢⱼ[1,1:3]

    for i = 1:3
        for j = 1:3
            SE.M.gᵢⱼ[j,i] += SE.M.Rᵢⱼ[1,j]*SE.M.Rᵢⱼ[1,i]
            for k = 1:3
                SE.M.dgᵢⱼ[i,j,k] += SE.M.Rᵢⱼ[k,j]*SE.M.Rᵢⱼ[1,i] + SE.M.Rᵢⱼ[1,j]*SE.M.Rᵢⱼ[k,i]
            end
        end
    end

    SE.M.gᵢⱼ[2,2] += T(1)
    SE.M.gᵢⱼ[3,3] += T(1)

    return 
end

"""
    R_Derivative!
Compute R and its derivatives in place
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

    if !P.StellaratorSymmetry
        t[3,:] = alss*Coord.Rbs[:,vᵢ] + blss*Coord.Rbs[:,vᵢ+1]
        t[4,:] = -0.5*Coord.Rbs[:,vᵢ] + 0.5*Coord.Rbs[:,vᵢ+1]
    end


    M.R₀₀       = dot(t[1,:],cosα)                      #R coordinate
    
    M.Rᵢⱼ[1,1]  = dot(t[2,:],cosα)                      #dRds
    M.Rᵢⱼ[1,2]  = dot(t[1,:],-Coord.mᵢ.*sinα)           #dRdu
    M.Rᵢⱼ[1,3]  = dot(t[1,:],Coord.nᵢ.*sinα)            #dRdv
    
    M.Rᵢⱼ[2,1]  = T(0)                                  #d2Rds2
    M.Rᵢⱼ[2,2]  = dot(t[2,:],-Coord.mᵢ.*sinα)           #d2Rdsdu
    M.Rᵢⱼ[2,3]  = dot(t[2,:],Coord.nᵢ.*sinα)            #d2Rdsdv
    
    M.Rᵢⱼ[3,2]  = dot(t[1,:],-Coord.mᵢ.^2 .*cosα)       #d2Rdu2
    M.Rᵢⱼ[3,3]  = dot(t[1,:],Coord.mᵢ.*Coord.nᵢ .*cosα) #d2Rdudv

    M.Rᵢⱼ[4,3]  = dot(t[1,:],-Coord.nᵢ.^2 .*cosα)       #d2Rdv2

    if !P.StellaratorSymmetry
        M.R₀₀       += dot(t[3,:],sinα)

        M.Rᵢⱼ[1,1]  += dot(t[4,:],sinα)
        M.Rᵢⱼ[1,2]  += dot(t[3,:],Coord.mᵢ.*cosα)
        M.Rᵢⱼ[1,3]  += dot(t[3,:],-Coord.nᵢ.*cosα)

        M.Rᵢⱼ[2,1]  += dot(t[3,:],sinα)
        M.Rᵢⱼ[2,2]  += dot(t[3,:],sinα)
        M.Rᵢⱼ[2,3]  += dot(t[3,:],sinα)

        M.Rᵢⱼ[3,2]  += dot(t[3,:],sinα)
        M.Rᵢⱼ[3,3]  += dot(t[3,:],sinα)
        M.Rᵢⱼ[4,3]  += dot(t[3,:],sinα)
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


