#linear basis is not high enough order for SAW, but useful for verification of simpler cases,
#such as the Helmholtz test case.

struct LB1d
    H :: Array{Float64, 2}
    dH :: Array{Float64, 2}
end

struct LB2d
    H :: Array{Float64, 4}
    dHx :: Array{Float64, 4}
    dHy :: Array{Float64, 4}
end


function linear_basis(gp::Array{Float64})

    S = LB1d(zeros(2, length(gp)), zeros(2, length(gp)))
    #don't think this is needed for linear.
    #t = @. (gp + 1)/(2) #converts to correct range for spline
    #H = zeros(2, length(gp))
    #dH = zeros(2, length(gp))

    S.H[1, :] = @. 1/2 * (1-gp)
    S.H[2, :] = @. 1/2 * (1+gp)

    S.dH[1, :] =  -1/2*ones(length(gp))
    S.dH[2, :] = 1/2*ones(length(gp))

    return S
end

function combine_basis(Sx1::LB1d, Sx2::LB1d, x1gp::Array{Float64}, x2gp::Array{Float64})

    as = (2, 2, length(x1gp), length(x2gp))

    S = LB2d(zeros(as), zeros(as), zeros(as))

    for i in 1:length(x1gp), j in 1:length(x2gp)

        for y in 1:2, x in 1:2

            S.H[x, y, i, j] = Sx1.H[x, i] * Sx2.H[y, j]
            S.dHx[x, y, i, j] = Sx1.dH[x, i] * Sx2.H[y, j]
            S.dHy[x, y, i, j] = Sx1.H[x, i] * Sx2.dH[y, j]
        end
    end

    return S

end

function linear_basis(x1gp::Array{Float64}, x2gp::Array{Float64})

    Sx = linear_basis(x1gp)
    Sy = linear_basis(x2gp)

    S = combine_basis(Sx, Sy, x1gp, x2gp)
    return S
end


function linear_local_to_global(xnode, ynode, ξx, ξy, xgrid, ygrid)
    
    if xnode == length(xgrid)
        dx = 2π + xgrid[1] - xgrid[xnode]
    else
        dx = xgrid[xnode+1] - xgrid[xnode]
    end

    if ynode == length(ygrid)
        dy = 2π + ygrid[1] - ygrid[ynode]
    else
        dy = ygrid[ynode+1] - ygrid[ynode]
    end

    mpx = @. (ξx + 1) /2 * dx
    mpy = @. (ξy + 1) /2 * dy

    xglobal = mpx .+ xgrid[xnode]
    yglobal = mpy .+ ygrid[ynode]

    return xglobal, yglobal, dx, dy
end

function linear_local_to_global(xnode, ynode, znode, ξx, ξy, ξz, xgrid, ygrid, zgrid)
    
    if xnode == length(xgrid)
        dx = 2π + xgrid[1] - xgrid[xnode]
    else
        dx = xgrid[xnode+1] - xgrid[xnode]
    end

    if ynode == length(ygrid)
        dy = 2π + ygrid[1] - ygrid[ynode]
    else
        dy = ygrid[ynode+1] - ygrid[ynode]
    end

    if znode == length(zgrid)
        dz = 2π + zgrid[1] - zgrid[znode]
    else
        dz = zgrid[znode+1] - zgrid[znode]
    end

    mpx = @. (ξx + 1) /2 * dx
    mpy = @. (ξy + 1) /2 * dy
    mpz = @. (ξz + 1) /2 * dz

    xglobal = mpx .+ xgrid[xnode]
    yglobal = mpy .+ ygrid[ynode]
    zglobal = mpz .+ zgrid[znode]

    return xglobal, yglobal, zglobal, dx, dy, dz
end

function linear_local_to_global(xnode, ξ, xgrid)
    
    dx = xgrid[xnode+1] - xgrid[xnode]
    
    mpx = @. (ξ + 1) /2 * dx

    xglobal = mpx .+ xgrid[xnode]

    return xglobal, dx
end
