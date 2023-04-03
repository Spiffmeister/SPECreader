

"""
    get_spec_field

Get the magnetic vector potential and/or magnetic field
"""
function get_spec_field(VecPot::VectorPotential,P::Parameters, s::T, θ::T, ξ::T, Mpol::Integer) where T

    A = zeros(T,3)
    gB = zeros(T,3)
    dgB = zeros(T,(3,3))


    if VecPot.Isingular
        s̄ = max((s+1.0)/2.0, 0.0)
        # Zernike()
    else
        # x = chebpoints()
        # c = chebinterp()
        poly = Chebychev(s,VecPot.Lrad)
    end


    for i = 1:VecPot.mn

        αᵢ = VecPot.mᵢ[i]*θ - VecPot.nᵢ[i]
        cosα = cos(αᵢ)
        sinα = sin(αᵢ)

        for l = 1:(VecPot.Lrad+1)
            
            if VecPot.isingular
            else
            end

            A[2] += VecPot.Ate[l,i] * poly[l,1]*cosα
            A[3] += VecPot.Ate[l,i] * poly[l,1]*cosα

            gB[1] += -(P.mᵢ[i] * VecPot.Aze[l,i] + P.nᵢ[i]*VecPot.Ate[l,i]) * poly[l,1] * sinα
            gB[2] += -VecPot.Aze[l,i] * poly[l,2] * cosα
            gB[3] += VecPot.Ate[l,i] * poly[l,2] * cosα

            ### DERIVATIVES ###
            # ds
            dgB[1,1] += -(P.mᵢ[i] * VecPot.Aze[l,i] + P.nᵢ[i] * VecPot.Ate[l,i]) * poly[l,2] * sinα
            dgB[2,1] += -VecPot.Aze[l,i] * poly[l,3] * cosα
            dgB[3,1] += VecPot.Ate[l,i] * poly[l,3] * cosα
            # dθ
            dgB[1,2] += -P.mᵢ[i] * (P.mᵢ[i] * VecPot.Aze[l,i] + P.nᵢ[i] * VecPot.Ate[l,i]) * poly[l,1] * cosα
            dgB[2,2] += P.mᵢ[i] * VecPot.Aze[l,i] * poly[l,2] * sinα
            dgB[3,2] += -P.mᵢ[i] * VecPot.Ate[l,i] * poly[l,2] * sinα
            # dζ
            dgB[1,1] += P.nᵢ[i] * (P.mᵢ[i] * VecPot.Aze[l,i] + P.nᵢ[i] * VecPot.Ate[l,i]) * poly[l,1] * cosα
            dgB[2,1] += -P.nᵢ[i] * VecPot.Aze[l,i] * poly[l,2] * sinα
            dgB[3,1] += P.nᵢ[i] * VecPot.Ate[l,i] * poly[l,2] * sinα


            if !P.StellaratorSymmetry
                A[2] += VecPot.Ato[l,i] * poly[l,1]*sinα
                A[3] += VecPot.Azo[l,i] * poly[l,1]*sinα
    
                gB[1] += -(P.mᵢ[i] * VecPot.Azo[l,i] + P.nᵢ[i]*VecPot.Ato[l,i]) * poly[l,1] * cosα
                gB[2] += -VecPot.Azo[l,i] * poly[l,2] * sinα
                gB[3] += VecPot.Ato[l,i] * poly[l,2] * sinα
    
                ### DERIVATIVES ###
                # ds
                dgB[1,1] += (P.mᵢ[i]*VecPot.Azo[l,i] + P.nᵢ[i]*VecPot.Ato[l,i]) * poly[l,2] * cosα
                dgB[2,1] += -VecPot.Azo[l,i] * poly[l,3] * sinα
                dgB[3,1] += VecPot.Ato[l,i] * poly[l,3] * sinα
                # dθ
                dgB[1,2] += -P.mᵢ[i] * (P.mᵢ[i]*VecPot.Azo[l,i] + P.nᵢ[i]*VecPot.Ato[l,i]) * poly[l,1] * sinα
                dgB[2,2] += -P.mᵢ[i] * VecPot.Aze[l,i] * poly[l,2] * cosα
                dgB[3,2] += P.mᵢ[i] * VecPot.Ate[l,i] * poly[l,2] * cosα
                # dζ
                dgB[1,1] += P.nᵢ[i] * (P.mᵢ[i]*VecPot.Azo[l,i] + P.nᵢ[i]*VecPot.Ato[l,i]) * poly[l,1] * sinα
                dgB[2,1] += P.nᵢ[i] * VecPot.Azo[l,i] * poly[l,2] * cosα
                dgB[3,1] += -P.nᵢ[i] * VecPot.Ato[l,i] * poly[l,2] * cosα
            end

        end

    end

end

