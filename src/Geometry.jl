
# """
#     GetMetric
# Get the metric information at (s,θ,ξ)
# """
# function GetMetric end
# function GetMetric(s::TT, θ::TT, ξ::TT, SpecVol::SPECequilibrium{TT}, vᵢ::Int, side::Symbol) where TT

#     if side == :inner
#         r_sign = -1/2
#         that_I = TT(lvol + 1)
#     else
#         r_sign = 1/2
#         that_I = TT(lvol)
#     end


#     α = mᵢ * θ - nᵢ * ζ
#     sinterm, costerm = sincos(α)

#     dR[1,1,1] = Rbc * costerm
#     dR[2,1,1] = Rbc * sinterm * -mᵢ
#     dR[3,1,1] = Rbc * sinterm * nᵢ

#     if SpecVol.CoordinateSingularity
#         dR[1,1,1] = Rbc * mᵢ * costerm + 2 * (Rbc - Rbc)*costerm/2
#     else
#         dR[1,1,1] = dR[1,1,1] - sum(Rbc * costerm) * 
#     end

# end






"""
    get_boundary(SpecVol::SPECEquilibrium{TT,ATT,TSA},lvol::Int) where {TT,ATT,TSA}

Get the boundary Fourier modes and return a named tuple with `Rbc`, `Zbs`, `m` and `n`.
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

    # BoundaryOut = Dict("Rbc" => Rbc, "Zbs" => Zbs, "m" => m, "n" => n)
    boundary = (Rbc=SpecVol.Rbc[:, 2],
            Zbs=SpecVol.Zbs[:, 2],
            m=SpecVol.m,
            n=SpecVol.n)
    
    return boundary
end


"""
    get_axis(SpecVol)

Get the Fourier modes of the magnetic axis of a SPEC equilibrium and output a named tuple `axis.R` and `axis.Z`.
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
                if SpecVol.Geometry == Cylindrical
                    fj = s̄
                else
                    fj = s̄^2
                end
            else
                if SpecVol.Geometry == Cylindrical
                    fj = s̄^(m + 1)
                else
                    fj = s̄^(m)
                end
            end


            Remn = Rbc[i, 1] + (Rbc[i, 2] - Rbc[i, 1]) * fj

            if !SpecVol.StellaratorSymmetric
                Romn = Rbs[i, 1] + (Rbs[i, 2] - Rbs[i, 1]) * fj
            end

            if SpecVol.Geometry == Toroidal
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
                Romn = zero(TT)
            end

            if SpecVol.Geometry == Toroidal
                Zomn = alss * Zbs[i, lvol] + blss * Zbs[i, lvol+1]
                if !SpecVol.StellaratorSymmetric
                    Zemn = alss * Zbc[i, lvol] + blss * Zbc[i, lvol+1]
                else
                    Zemn = zero(TT)
                end
            else
                Zomn = zero(TT)
                Zemn = zero(TT)
            end

        end

        R += Remn * costerm + Romn * sinterm
        if SpecVol.Geometry == Toroidal
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
    tmp = get_RZ(q[1], q[2], ζ, SpecVol, lvol)
    tmp = norm(X .- tmp)
    return tmp
end
function _coordinate_difference!(qout, q, X, ζ, SpecVol::SPECEquilibrium, lvol::Integer)
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
