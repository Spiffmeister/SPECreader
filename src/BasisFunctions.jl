

"""
    get_zernike!(zernike, s::TT, SpecVol::SPECEquilibrium) where TT

Get the Zernike weights at the point ``s\\in[-1,1]`` in a given SPEC equilibrium.
"""
function get_zernike!(zernike, s::TT, SpecVol::SPECEquilibrium,derivative=false) where {TT}

    RadialResolution = SpecVol.RadialResolution[1]

    zernike .= zero(TT)

    PoloidalResolution = SpecVol.PoloidalResolution


    s̄ = (1 + s) / 2 # in SPEC this is s̄

    rᵐ = one(TT)
    rᵐ⁻¹ = zero(TT)
    rᵐ⁻² = zero(TT)

    for m in 0:PoloidalResolution

        mᵢ = m + 1 # Adjusted for indexing

        if RadialResolution ≥ m
            zernike[mᵢ, mᵢ, 1] = rᵐ
            zernike[mᵢ, mᵢ, 2] = TT(m) * rᵐ⁻¹
            if derivative
                zernike[mᵢ,mᵢ,3] = TT(m * (m-1)) * rᵐ⁻²
            end
        end

        if RadialResolution ≥ m + 2
            zernike[mᵢ+2, mᵢ, 1] = TT(m + 2) * rᵐ * s̄^2 - TT(m + 1) * rᵐ
            zernike[mᵢ+2, mᵢ, 2] = TT(m + 2)^2 * rᵐ * s̄ - TT(m + 1) * TT(m) * rᵐ⁻¹
            if derivative
                zernike[mᵢ+2,mᵢ,3] = TT(m+2)^2 * TT(m+1) * rᵐ⁻¹ - TT(m+1) * TT(m) * TT(m-1) * rᵐ⁻²
            end
        end

        for n in m+4:2:RadialResolution
            nᵢ = n + 1
            factor1 = TT(n) / TT(n^2 - m^2)
            factor2 = TT(4 * (n - 1))
            factor3 = TT(n - 2 + m)^2 / TT(n - 2) + TT(n - m)^2 / TT(n)
            factor4 = TT((n - 2)^2 - m^2) / TT(n - 2)

            zernike[nᵢ, mᵢ, 1] = factor1 * ((factor2 * s̄^2 - factor3) * zernike[nᵢ-2, mᵢ, 1] - factor4 * zernike[nᵢ-4, mᵢ, 1])
            zernike[nᵢ, mᵢ, 2] = factor1 * (2 * factor2 * s̄ * zernike[nᵢ-2, mᵢ, 1] + (factor2 * s̄^2 - factor3) * zernike[nᵢ-2, mᵢ, 2] - factor4 * zernike[nᵢ-4, mᵢ, 2])
            if derivative
                zernike[nᵢ,mᵢ,3] = factor1 * (2 * factor2 * (2 * rᵐ * zernike[nᵢ-2,mᵢ,1] + zernike[nᵢ-2,mᵢ,1])) + 
                    (facto2 * rᵐ^2 - factor3) * zernike[nᵢ-2,mᵢ,3] - factor4*zernike[nᵢ-4,mᵢ,3]
            end
        end

        rᵐ⁻¹ = rᵐ
        rᵐ = rᵐ * s̄

    end

    for n in 2:2:RadialResolution
        nᵢ = n + 1
        zernike[nᵢ, 1, 1] = zernike[nᵢ, 1, 1] - (-1)^(TT(n) / 2)
    end


    if PoloidalResolution ≥ 1
        for n in 3:2:RadialResolution
            nᵢ = n + 1
            zernike[nᵢ, 2, 1] = zernike[nᵢ, 2, 1] - (-1)^(TT(n - 1) / 2) * TT(n + 1) / 2 * s̄
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
    get_chebychev!(chebychev,s::TT,SpecVol::SPECEquilibrium,lvol::Int) where TT

Chebyshev coefficients are provided by Rbc, Rbs, Zbc, Zbs
"""
function get_chebychev!(chebychev,s::TT,SpecVol::SPECEquilibrium,lvol::Int,derivative=false) where TT
    
    RadialResolution = SpecVol.RadialResolution[lvol]

    chebychev .= zero(TT)

    chebychev[1,1] = one(TT)
    # chebychev[1,2] = zero(TT)
    chebychev[2,1] = s
    chebychev[2,2] = one(TT)

    if derivative
        # chebychev[1,3] = zero(TT)
        # chebychev[2,3] = zero(TT)
    end

    for l in 3:RadialResolution+1
        chebychev[l,1] = 2 * s * chebychev[l-1,1] - chebychev[l-2,1]
        chebychev[l,2] = 2 * chebychev[l-1,1] + 2 * s * chebychev[l-1,2] - chebychev[l-2,2]
        if derivative
            chebychev[l,3] = 2 * chebychev[l-1,2] + 2 * chebychev[l-1,2] + 2 * s * chebychev[l-1,3] - chebychev[l-2,3]
        end
    end

    for l in 2:RadialResolution+1
        chebychev[l,1] = chebychev[l,1] - (-1)^TT(l-1)
    end

    for l in 1:RadialResolution+1
        chebychev[l,1] = chebychev[l,1] / TT(l)
        chebychev[l,2] = chebychev[l,2] / TT(l)
    end

end
