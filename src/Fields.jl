


"""
    get_Bfield(s::TT, θ::TT, ζ::TT, SpecVol::SPECEquilibrium{TT,ATT,TSA},lvol::Int=1) where {TT,ATT,TSA}

Take a point in ``(s,\\theta,\\zeta)`` and return the covarient (cotangent) components of the magnetic field.
"""
function get_Bfield(s::TT, θ::TT, ζ::TT, SpecVol::SPECEquilibrium{TT,ATT,TSA},lvol::Int=1) where {TT,ATT,TSA}

    # B = SpecVol.cache
    # B .= zero(TT)

    B1 = zero(TT)
    B2 = zero(TT)
    B3 = zero(TT)


    mn = SpecVol.mn
    m = SpecVol.m
    n = SpecVol.n

    Ate = SpecVol.Ate
    Aze = SpecVol.Aze
    Ato = SpecVol.Ato
    Azo = SpecVol.Azo

    RadialResolution = SpecVol.RadialResolution[lvol]

    if lvol == 1
        CoordinateSingularity = SpecVol.CoordinateSingularity
        RadialOffset = 0
    else
        CoordinateSingularity = false
        RadialOffset = sum(@views SpecVol.RadialResolution[1:lvol-1]) + 1
    end

    RadialBasis = SpecVol.RadialBasis::TSA
    if CoordinateSingularity
        get_zernike!(RadialBasis, s, SpecVol)
    else
        chebychev = view(RadialBasis, 1:RadialResolution+1, 1:2, 1)
        get_chebychev!(chebychev, s, SpecVol, lvol)
    end

    for i = 1:mn

        mᵢ = TT(m[i])
        nᵢ = TT(n[i])

        mindex = m[i] + 1 # m starts from 0

        if CoordinateSingularity
            Poly1 = view(RadialBasis, 1:RadialResolution+1, mindex, 1)
            Poly2 = view(RadialBasis, 1:RadialResolution+1, mindex, 2)
        else
            Poly1 = view(RadialBasis, 1:RadialResolution+1, 1, 1)
            Poly2 = view(RadialBasis, 1:RadialResolution+1, 2, 1)
        end

        α = mᵢ * θ - nᵢ * ζ
        sterm, cterm = sincos(α)

        for l in 1:RadialResolution+1
            ll = RadialOffset + l
            B1 += (-mᵢ * Aze[ll, i] - nᵢ * Ate[ll, i]) * Poly1[l] * sterm
            B2 += -Aze[ll, i] * Poly2[l] * cterm
            B3 += Ate[ll, i] * Poly2[l] * cterm
        end

        if SpecVol.StellaratorSymmetric
            for l in 1:RadialResolution+1
                ll = RadialOffset + l
                B1 += (mᵢ * Azo[ll, i] + nᵢ * Ato[ll, i]) * Poly1[l] * cterm
                B2 += -Azo[ll, i] * Poly2[l] * sterm
                B3 += Ato[ll, i] * Poly2[l] * sterm
            end
        end
    end

    return B1, B2, B3
end




"""
    field_line!(ẋ, t, x, SpecVol::SPECEquilibrium,lvol::Int=1)

In place field line tracing function given a `lvol`

TODO: Implement for multi-volume spec equilibria
"""
function field_line!(ẋ, t, x, SpecVol::SPECEquilibrium,lvol::Int=1)

    B1, B2, B3 = get_Bfield(x[1], x[2], t, SpecVol,lvol)

    ẋ[1] = B1 / B3
    ẋ[2] = B2 / B3

end
