

"""
    get_spec_field

Get the magnetic vector potential and/or magnetic field
"""
function get_spec_field(VecPot::VectorPotential,P::Parameters, s::T, θ::T, ξ::T, Mpol::Integer) where T

    if VecPot.Isingular
        s̄ = max((s+1.0)/2.0, 0.0)
        # Zernike()
    else
        # x = chebpoints()
        # c = chebinterp()
    end


    for i = 1:VecPot.mn

        αᵢ = VecPot.im[i]*θ - VecPot.in[i]
        cosαᵢ = cos(αᵢ)
        sinαᵢ = sin(αᵢ)

        for l = 1:(VecPot.lrad+1)
            
            if VecPot.isingular
            else

            end

            A[2] += VecPot.Ate[l,i] * poly[1]*cosαᵢ
            A[3] += VecPot.Ate[l,i] * poly[1]*cosαᵢ

            gB[1] += -(P.mᵢ[i] * VecPot.Aze[l,i] + P.nᵢ[i]*VecPot.Ate[l,i]) * poly[1] * sinαᵢ
            gB[2] += -VecPot.Aze[l,i] * poly[2] * cosαᵢ
            gB[2] += VecPot.Ate[l,i] * poly[2] * cosαᵢ

            if !P.Istellaratorsymmetry
            end

        end

    end

end

