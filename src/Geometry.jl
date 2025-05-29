
"""
    GetMetric
Get the metric information at (s,θ,ξ)
"""
function GetMetric end
function GetMetric(SE::SPECequilibrium{Cartesian,V}, s::T, θ::T, ξ::T, vᵢ::Integer) where {T,V}

    R_Derivative!(SE.M, s, θ, ξ, vᵢ, SE.Ri, SE.Params)

    SE.M.x[1] = θ
    SE.M.x[2] = ξ
    SE.M.x[3] = SE.M.R₀₀

    SE.M.gᵢⱼ[:, :] .= T(0)
    SE.M.dgᵢⱼ[:, :, :] .= T(0)

    SE.M.Jac = SE.M.Rᵢⱼ[1, 1]
    SE.M.∇Jac = SE.M.Rᵢⱼ[1:3, 1]

    SE.M.JacMat[1, 2] = T(1)
    SE.M.JacMat[2, 3] = T(1)
    SE.M.JacMat[3, 1:3] = SE.M.Rᵢⱼ[1, 1:3]

    for i = 1:3
        for j = 1:3
            SE.M.gᵢⱼ[j, i] += SE.M.Rᵢⱼ[1, j] * SE.M.Rᵢⱼ[1, i]
            for k = 1:3
                SE.M.dgᵢⱼ[i, j, k] += SE.M.Rᵢⱼ[k, j] * SE.M.Rᵢⱼ[1, i] + SE.M.Rᵢⱼ[1, j] * SE.M.Rᵢⱼ[k, i]
            end
        end
    end

    SE.M.gᵢⱼ[2, 2] += T(1)
    SE.M.gᵢⱼ[3, 3] += T(1)

    return
end

"""
    R_Derivative!
Compute R and its derivatives in place
"""
function R_Derivative!(M::Metrics, s::T, θ::T, ξ::T, vᵢ::Integer, Coord::Coordinates, P::Parameters) where {T}

    a = αᵢ(Coord.mᵢ, Coord.nᵢ, θ, ξ)
    cosα = cos.(a)
    sinα = sin.(a)

    t = zeros(T, (4, P.mn))

    alss = (1.0 - s) / 2.0
    blss = (1.0 + s) / 2.0

    t[1, :] = alss * Coord.Rbc[:, vᵢ] + blss * Coord.Rbc[:, vᵢ+1]
    t[2, :] = -0.5 * Coord.Rbc[:, vᵢ] + 0.5 * Coord.Rbc[:, vᵢ+1]

    if !P.StellaratorSymmetry
        t[3, :] = alss * Coord.Rbs[:, vᵢ] + blss * Coord.Rbs[:, vᵢ+1]
        t[4, :] = -0.5 * Coord.Rbs[:, vᵢ] + 0.5 * Coord.Rbs[:, vᵢ+1]
    end


    M.R₀₀ = dot(t[1, :], cosα)                      #R coordinate

    M.Rᵢⱼ[1, 1] = dot(t[2, :], cosα)                      #dRds
    M.Rᵢⱼ[1, 2] = dot(t[1, :], -Coord.mᵢ .* sinα)           #dRdu
    M.Rᵢⱼ[1, 3] = dot(t[1, :], Coord.nᵢ .* sinα)            #dRdv

    M.Rᵢⱼ[2, 1] = T(0)                                  #d2Rds2
    M.Rᵢⱼ[2, 2] = dot(t[2, :], -Coord.mᵢ .* sinα)           #d2Rdsdu
    M.Rᵢⱼ[2, 3] = dot(t[2, :], Coord.nᵢ .* sinα)            #d2Rdsdv

    M.Rᵢⱼ[3, 2] = dot(t[1, :], -Coord.mᵢ .^ 2 .* cosα)       #d2Rdu2
    M.Rᵢⱼ[3, 3] = dot(t[1, :], Coord.mᵢ .* Coord.nᵢ .* cosα) #d2Rdudv

    M.Rᵢⱼ[4, 3] = dot(t[1, :], -Coord.nᵢ .^ 2 .* cosα)       #d2Rdv2

    if !P.StellaratorSymmetry
        M.R₀₀ += dot(t[3, :], sinα)

        M.Rᵢⱼ[1, 1] += dot(t[4, :], sinα)
        M.Rᵢⱼ[1, 2] += dot(t[3, :], Coord.mᵢ .* cosα)
        M.Rᵢⱼ[1, 3] += dot(t[3, :], -Coord.nᵢ .* cosα)

        M.Rᵢⱼ[2, 1] += dot(t[3, :], sinα)
        M.Rᵢⱼ[2, 2] += dot(t[3, :], sinα)
        M.Rᵢⱼ[2, 3] += dot(t[3, :], sinα)

        M.Rᵢⱼ[3, 2] += dot(t[3, :], sinα)
        M.Rᵢⱼ[3, 3] += dot(t[3, :], sinα)
        M.Rᵢⱼ[4, 3] += dot(t[3, :], sinα)
    end


    # Rᵢⱼᵀ = Rᵢⱼ
    M.Rᵢⱼ[3, 1] = M.Rᵢⱼ[2, 2]
    M.Rᵢⱼ[4, 1] = M.Rᵢⱼ[2, 3]
    M.Rᵢⱼ[4, 2] = M.Rᵢⱼ[3, 3]

end

"""
    Z_Derivative
"""
function Z_Derivative()

    # cosa = cos.(αᵢ())
    # sina = sin.(αᵢ())

end


@inline αᵢ(m, n, θ, ξ) = m * θ - n * ξ





"""
    get_boundary(SpecVol::SPECEquilibrium{TT,ATT,TSA},lvol::Int) where {TT,ATT,TSA}

Get the boundary.
"""
function get_boundary(SpecVol::SPECEquilibrium{TT,ATT,TSA}, lvol::Int) where {TT,ATT,TSA}
    if lvol == 1
        boundary = (Rbc=SpecVol.Rbc[:, 2],
            Zbs=SpecVol.Zbs[:, 2],
            m=SpecVol.m,
            n=SpecVol.n)
    else
        error("lvol>1 not implemented yet")
    end
    return boundary
end
"""
    get_boundary(SpecVol::SPECEquilibrium;lvol=SpecVol.NumberofVolumes)

Pull the outer boundary from a spec equlibrium data structure.

TODO: Add functionality to `lvol` so that we can pull the interface positions.
"""
function get_boundary(SpecVol::SPECEquilibrium; lvol=SpecVol.NumberofVolumes)
    Rbc = SpecVol.Rbc[:, 2]
    Zbs = SpecVol.Zbs[:, 2]

    PoloidalResolution = SpecVol.PoloidalResolution
    Ntor = SpecVol.Ntor

    m = vcat([collect(-PoloidalResolution:PoloidalResolution) for i in 1:Ntor+1]...)[PoloidalResolution+1:end]
    n = vcat(zeros(Int, PoloidalResolution + 1), repeat(1:Ntor, inner=2PoloidalResolution + 1))

    if eltype(m) == Int32
        m = [promote(m..., Int64(1))[1:end-1]...]
    end

    BoundaryOut = Dict("Rbc" => Rbc, "Zbs" => Zbs, "m" => m, "n" => n)
    return BoundaryOut
end


"""
    get_axis(SpecVol)

Get the magnetic axis.
"""
function get_axis(SpecVol)
    return (R=SpecVol.Rbc[:, 1], Z=SpecVol.Zbs[:, 1])
end





"""
    get_RZ(s::TT,θ,ζ,SpecVol::SPECEquilibrium,lvol::Integer) where TT

Get the ``(R,Z)`` coordinates from ``(s,\\theta,\\zeta)\\in[-1,1]\\times[0,2\\pi)\\times[0,2\\pi)`` logical coordinates.
"""
function get_RZ(s::TT, θ, ζ, SpecVol::SPECEquilibrium, lvol::Integer) where {TT}


    Rbc = SpecVol.Rbc
    Rbs = SpecVol.Rbs
    Zbc = SpecVol.Zbc
    Zbs = SpecVol.Zbs

    s̄ = (s + TT(1)) / TT(2)

    alss = 1 / 2 - s / 2 # TODO: terrible variable name
    blss = 1 / 2 + s / 2 # TODO: terrible variable name

    local R = zero(TT)
    local Z = zero(TT)

    for i in 1:SpecVol.mn

        m = SpecVol.m[i]
        n = SpecVol.n[i]

        sinterm, costerm = sincos(m * θ - n * ζ)

        Romn = zero(TT)
        Remn = zero(TT)
        Zomn = zero(TT)
        Zemn = zero(TT)

        if SpecVol.CoordinateSingularity

            if m == 0
                if SpecVol.Geometry == :cylindrical
                    fj = s̄
                else
                    fj = s̄^2
                end
            else
                if SpecVol.Geometry == :cylindrical
                    fj = s̄^(m + 1)
                else
                    fj = s̄^(m)
                end
            end


            Remn = Rbc[i, 1] + (Rbc[i, 2] - Rbc[i, 1]) * fj

            if !SpecVol.StellaratorSymmetric
                Romn = Rbs[i, 1] + (Rbs[i, 2] - Rbs[i, 1]) * fj
            end

            if SpecVol.Geometry == :toroidal
                Zomn = Zbs[i, 1] + (Zbs[i, 2] - Zbs[i, 1]) * fj
                if !SpecVol.StellaratorSymmetric
                    Zemn = Zbc[i, 1] + (Zbc[i, 2] - Zbc[i, 1]) * fj
                end
            end

        else

            Remn = alss * Rbc[i, lvol] + blss * Rbc[i, lvol+1]

            if !SpecVol.StellaratorSymmetric
                Romn = alss * Rbs[i, lvol] + blss * Rbs[i, lvol+1]
            else
                Romn = TT(0)
            end

            if SpecVol.Geometry == :toroidal
                Zomn = alss * Zbs[i, lvol] + blss * Zbs[i, lvol+1]
                if !SpecVol.StellaratorSymmetric
                    Zemn = alss * Zbc[i, lvol] + blss * Zbc[i, lvol+1]
                else
                    Zemn = TT(0)
                end
            else
                Zomn = TT(0)
                Zemn = TT(0)
            end

        end

        R += Remn * costerm + Romn * sinterm
        if SpecVol.Geometry == :toroidal
            Z += Zemn * costerm + Zomn * sinterm
        end

    end

    return R, Z

end




"""
    _coordinate_difference(q,X,ζ,SpecVol::SPECEquilibrium,lvol::Integer)

Finding ``(s,\\theta,\\zeta)`` points for field line tracing. ``X = (R,Z)`` and ``q = (s,\\theta)``, ``\\zeta`` is the toroidal angle
"""
function _coordinate_difference(q, X, ζ, SpecVol::SPECEquilibrium, lvol::Integer)
    # X
    tmp = get_RZ(q[1], q[2], ζ, SpecVol, lvol)
    tmp = norm(X .- tmp)
    return tmp
end
function _coordinate_difference!(qout, q, X, ζ, SpecVol::SPECEquilibrium, lvol::Integer)
    # X
    qout .= get_RZ(q[1], q[2], ζ, SpecVol, lvol)
    qout .= X .- qout
end



"""
    find_sθζ(X,ζ,SpecVol::SPECEquilibrium,lvol::Integer,max_attempts=100)

Find the ``(s,\\theta,\\zeta)`` point corresponding to a given ``(R,Z)`` point for a fixed ``\\zeta``.
It is possible that the starting location is bad and the located point is outside of the computational domain. If we generate random initial conditions until it ends up inside the domain.
"""
function find_sθζ(X, ζ, SpecVol::SPECEquilibrium, lvol::Integer, max_attempts=100, tol=1e-12)

    sol = optimize(q -> _coordinate_difference(q, X, ζ, SpecVol, lvol), [0.0, π])
    u = sol.minimizer

    # This ensures things are inside the box, probably need to check that it doesn't cause some weird issues with point initialisation
    attempt = 1
    while !(-1 ≤ u[1] ≤ 1) || !(0.0 ≤ u[2] ≤ 2π) || (attempt > max_attempts)
        IC = rand(2)
        IC[1] = 2IC[2] - 1
        IC[2] = IC[2] * 2π
        sol = optimize(q -> _coordinate_difference(q, X, ζ, SpecVol, lvol), IC)
        u = sol.minimizer
        attempt += 1
    end

    nlp = NonlinearProblem((qout, q, X) -> _coordinate_difference!(qout, q, X, ζ, SpecVol, lvol), u, X)
    u = solve(nlp).u

    return u
end
