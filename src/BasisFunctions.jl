





struct zernike{T}
    z::AbstractArray{T}
    Dz::AbstractArray{T}
    DDz::AbstractArray{T}
end



"""
    get_zernike!(zernike, s::TT, SpecVol::SPECEquilibrium) where TT

Get the Zernike weights at the point ``s\\in[-1,1]`` in a given SPEC equilibrium.
"""
function get_zernike!(zernike, s::TT, SpecVol::SPECEquilibrium) where {TT}

    RadialResolution = SpecVol.RadialResolution[1]

    zernike .= zero(TT)

    PoloidalResolution = SpecVol.PoloidalResolution


    s = (1 + s) / 2 # in SPEC this is s̄

    rᵐ = one(TT)
    rᵐ⁻¹ = zero(TT)

    for m in 0:PoloidalResolution

        mᵢ = m + 1 # Adjusted for indexing

        if RadialResolution ≥ m
            zernike[mᵢ, mᵢ, 1] = rᵐ
            zernike[mᵢ, mᵢ, 2] = TT(m) * rᵐ⁻¹
        end

        if RadialResolution ≥ m + 2
            zernike[mᵢ+2, mᵢ, 1] = TT(m + 2) * rᵐ * s^2 - TT(m + 1) * rᵐ
            zernike[mᵢ+2, mᵢ, 2] = TT(m + 2)^2 * rᵐ * s - TT(m + 1) * TT(m) * rᵐ⁻¹
        end

        for n in m+4:2:RadialResolution
            nᵢ = n + 1
            factor1 = TT(n) / TT(n^2 - m^2)
            factor2 = TT(4 * (n - 1))
            factor3 = TT(n - 2 + m)^2 / TT(n - 2) + TT(n - m)^2 / TT(n)
            factor4 = TT((n - 2)^2 - m^2) / TT(n - 2)

            zernike[nᵢ, mᵢ, 1] = factor1 * ((factor2 * s^2 - factor3) * zernike[nᵢ-2, mᵢ, 1] - factor4 * zernike[nᵢ-4, mᵢ, 1])
            zernike[nᵢ, mᵢ, 2] = factor1 * (2 * factor2 * s * zernike[nᵢ-2, mᵢ, 1] + (factor2 * s^2 - factor3) * zernike[nᵢ-2, mᵢ, 2] - factor4 * zernike[nᵢ-4, mᵢ, 2])
        end

        rᵐ⁻¹ = rᵐ
        rᵐ = rᵐ * s

    end

    for n in 2:2:RadialResolution
        nᵢ = n + 1
        zernike[nᵢ, 1, 1] = zernike[nᵢ, 1, 1] - (-1)^(TT(n) / 2)
    end


    if PoloidalResolution ≥ 1
        for n in 3:2:RadialResolution
            nᵢ = n + 1
            zernike[nᵢ, 2, 1] = zernike[nᵢ, 2, 1] - (-1)^(TT(n - 1) / 2) * TT(n + 1) / 2 * s
            zernike[nᵢ, 2, 2] = zernike[nᵢ, 2, 2] - (-1)^(TT(n - 1) / 2) * TT(n + 1) / 2
        end
    end


    for m in 0:PoloidalResolution
        mᵢ = m + 1
        for n in m:2:RadialResolution
            nᵢ = n + 1
            zernike[nᵢ, mᵢ, 1] = zernike[nᵢ, mᵢ, 1] / TT(n + 1)
            zernike[nᵢ, mᵢ, 2] = zernike[nᵢ, mᵢ, 2] / TT(n + 1)
        end
    end

    for j in 1:PoloidalResolution+1
        for i in 1:RadialResolution+1
            zernike[i, j, 2] /= 2
        end
    end

end



"""
    Cheby
Chebyshev coefficients are provided by Rbc, Rbs, Zbc, Zbs
"""


function Cheby()

end
