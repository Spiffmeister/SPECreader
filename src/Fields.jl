


"""
    get_Bfield(s::TT,θ::TT,ζ::TT,SpecVol::SPECEquilibrium{TT,ATT,TSA}) where {TT,ATT,TSA}

Take a point in ``(s,\\theta,\\zeta)`` and return the covarient (cotangent) components of the magnetic field.
"""
function get_Bfield(s::TT, θ::TT, ζ::TT, SpecVol::SPECEquilibrium{TT,ATT,TSA}) where {TT,ATT,TSA}

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

    RadialResolution = SpecVol.RadialResolution[1]
    CoordinateSingularity = SpecVol.CoordinateSingularity


    if CoordinateSingularity
        RadialBasis = SpecVol.RadialBasis::TSA

        get_zernike!(RadialBasis, s, SpecVol)
    else

        raise("Cheby not implemented")
    end

    for i = 1:mn

        mᵢ = TT(m[i])
        nᵢ = TT(n[i])

        mindex = m[i] + 1 # m starts from 0

        if CoordinateSingularity
            Poly1 = view(RadialBasis, 1:RadialResolution+1, mindex, 1)
            Poly2 = view(RadialBasis, 1:RadialResolution+1, mindex, 2)
        else
        end

        α = mᵢ * θ - nᵢ * ζ
        sterm, cterm = sincos(α)

        for l in 1:RadialResolution+1
            # B[1] += (-mᵢ * Aze[l, i] - nᵢ * Ate[l, i])  * Poly1[l] * sterm
            # B[2] += - Aze[l, i]                         * Poly2[l] * cterm
            # B[3] += Ate[l, i]                           * Poly2[l] * cterm

            B1 += (-mᵢ * Aze[l, i] - nᵢ * Ate[l, i]) * Poly1[l] * sterm
            B2 += -Aze[l, i] * Poly2[l] * cterm
            B3 += Ate[l, i] * Poly2[l] * cterm
        end

        if SpecVol.StellaratorSymmetric
            for l in 1:RadialResolution+1
                # B[1] += (mᵢ * Azo[l, i] + nᵢ * Ato[l, i])   * Poly1[l] * cterm
                # B[2] += - Azo[l, i]                         * Poly2[l] * sterm
                # B[3] += Ato[l, i]                           * Poly2[l] * sterm

                B1 += (mᵢ * Azo[l, i] + nᵢ * Ato[l, i]) * Poly1[l] * cterm
                B2 += -Azo[l, i] * Poly2[l] * sterm
                B3 += Ato[l, i] * Poly2[l] * sterm
            end
        end
    end

    return B1, B2, B3
end




"""
  field_line!(ẋ, t, x, SpecVol::SPECEquilibrium)

In place field line tracing function given a `SpecVol`

TODO: Implement for multi-volume spec equilibria
"""
function field_line!(ẋ, t, x, SpecVol::SPECEquilibrium)

    B1, B2, B3 = get_Bfield(x[1], x[2], t, SpecVol)

    ẋ[1] = B1 / B3
    ẋ[2] = B2 / B3

end
