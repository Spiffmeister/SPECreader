





struct zernike{T}
    z   :: AbstractArray{T}
    Dz  :: AbstractArray{T}
    DDz :: AbstractArray{T}
end


"""
    zernike

Get Zernike polynomial
"""
function build_zernike(Z::zernike,r,lrad,mpol)

    rm = 1.0
    rm1 = 0.0
    rm2 = 0.0

    for m = 1:mpol+1

        if lrad ≥ m
            Z.z[m,m]    = rm
            Z.Dz[m,m]   = m*rm1
            Z.DDz[m,m]  = m*(m-1)*rm2
        end
        if lrad ≥ m+2
            Z.z[m+2,m]  = (m+2)*rm^2 - (m+1)*rm
            Z.Dz[m+2,m] = ((m+2)^2)*rm*r - (m+1)*m*rm1
            Z.DDz[m+2,m]= ((m+2)^2)*(m+1)*rm - (m+1)*m*(m-1)*rm2
        end

        for n = m+4:2:lrad
            
        end

        rm2 = rm1
        rm1 = rm
        rm = rm * r
    end

end

"""
    Cheby

"""
function Chebychev(s::T,Lrad::Integer) where T
    Lrad += 1
    cheb = zeros(T,(Lrad,3))
    cheb[1,1] = cheb[2,2] = T(1)
    cheb[2,1] = s

    for i = 3:Lrad
        cheb[i,1] = T(2)*s*cheb[i-1,1] - cheb[i-2,1]
        cheb[i,2] = T(2)*cheb[i-1,1] + T(2)*s*cheb[i-1,2] - cheb(i-2,2)
        cheb[i,3] = T(4)*cheb[i-1,2] + T(2)*s*cheb[i-1,3] - cheb(i-2,3)
    end

    for i = 2:Lrad
        cheb[i,1] -= T(-1)^i
    end

    for i = 1:Lrad
        cheb[i,:] ./= T(i+1)
    end

end



