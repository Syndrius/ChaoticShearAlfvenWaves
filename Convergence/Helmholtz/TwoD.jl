
function cubic_mat_size(Nx::Int64, Ny::Int64)
    return 4*Nx*Ny
end

function compute_linear_bcs(Nx::Int64, Ny::Int64)
    leftx = 1:Ny
    rightx = (Nx-1)*Ny:Nx*Ny

    lefty = 1:Ny:Nx*Ny

    righty = Nx:Ny:Nx*Ny

    #[x1y1, x1y2, ... x1, yN, x2y1, x2y2, .... xNy1, xNy2, ..., xNyN]
    #return collect(rightx)

    return vcat(leftx, rightx)
    return vcat(leftx, lefty, rightx, righty)
end

function compute_cubic_bcs(Nx::Int64, Ny::Int64)
    left_bc1 = 1:4:4*Ny
    left_bc2 = 2:4:4*Ny

    right_bc1 = 1 + (Nx - 1)*4*Ny:4:4*Nx*Ny
    right_bc2 = 2 + (Nx - 1)*4*Ny:4:4*Nx*Ny

    #return vcat(right_bc1, right_bc2)

    x_bc = vcat(left_bc1, left_bc2, right_bc1, right_bc2)

    #return x_bc #for periodic case.

    left_bc1 = 1:4*Ny:4*Nx*Ny
    left_bc2 = 3:4*Ny:4*Nx*Ny

    right_bc1 = 1+4*(Ny-1):4*Ny:4*Ny*Nx
    right_bc2 = 3+4*(Ny-1):4*Ny:4*Ny*Nx #oof with this one!

    y_bc = vcat(left_bc1, left_bc2, right_bc1, right_bc2)


    #ooft bcs in y are going to be cooked af.

    return vcat(x_bc, y_bc)

end

function cubic_bcs_cart(Nx::Int64, Ny::Int64)

    bcs = Int64[]
    for i in 1:4*Nx*Ny
        xind, yind, h = index_to_grid_cubic(i, Nx, Ny)
        if h == 1 || h == 2
            if xind == 1 || xind == Nx
                push!(bcs, i)
            end
        end
        if h == 1 || h == 3
            if yind == 1 || yind == Ny
                push!(bcs, i)
            end
        end
    end
    return bcs

end
#start with cartesian case, we will want to shift to cylindrical later to more accuratly represent our real equation.
function linear(Nx, Ny, gpx, gpy)


    xgrid = LinRange(0, π, Nx)#+1)[1:end-1]
    ygrid = LinRange(0, π, Ny)#+1)[1:end-1]
    xgrid = LinRange(0, 1, Nx)#+1)[1:end-1]
    ygrid = LinRange(0, 2π, Ny+1)[1:end-1]

    ξx, wgx = gausslegendre(gpx)
    ξy, wgy = gausslegendre(gpy)

    S = linear_basis(ξx, ξy)

    M = zeros(Float64, Nx*Ny, Nx*Ny)

    A = zeros(Float64, Nx*Ny, Nx*Ny)

    bcs = compute_linear_bcs(Nx, Ny)
    

    #for i in 1:Nx-1, j in 1:Ny-1
    for i in 1:Nx-1, j in 1:Ny

        x, y, dx, dy = linear_local_to_global(i, j, ξx, ξy, xgrid, ygrid)

        jac = dx * dy / 4

        for trialx in 1:2, trialy in 1:2

            r_ind = grid_to_index_linear(i, j, trialx, trialy, Nx, Ny)

            if r_ind in bcs
                continue
            end

            for testx in 1:2, testy in 1:2

                l_ind = grid_to_index_linear(i, j, testx, testy, Nx, Ny)

                if l_ind in bcs
                    continue
                end

                for kx in 1:gpx, ky in 1:gpy

                    #M[l_ind, r_ind] += jac * (S.dHx[testx, testy, kx, ky] * S.dHx[trialx, trialy, kx, ky] * 2/dx * 2/dx + S.dHy[testx, testy, kx, ky] * S.dHy[trialx, trialy, kx, ky] * 2/dy * 2/dy) * wgx[kx] * wgy[ky]

                    #display(jac * (S.dHx[testx, testy, kx, ky] * S.dHx[trialx, trialy, kx, ky] + S.dHy[testx, testy, kx, ky] * S.dHy[trialx, trialy, kx, ky]) / jac^2 * wgx[kx] * wgy[ky])

                    #A[l_ind, r_ind] += jac * S.H[testx, testy, kx, ky] * S.H[trialx, trialy, kx, ky] * wgx[kx] * wgy[ky]

                    #polar Helmholtz
                    M[l_ind, r_ind] += jac * (+S.dHx[testx, testy, kx, ky] * S.dHx[trialx, trialy, kx, ky] * (2/dx)^2 
                                              - 1/x[kx] * S.H[testx, testy, kx, ky] * S.dHx[trialx, trialy, kx, ky] * (2/dx)
                                              + 1/x[kx]^2 * S.dHy[testx, testy, kx, ky] * S.dHy[trialx, trialy, kx, ky] * (2/dy)^2) * wgx[kx] * wgy[ky]

                    A[l_ind, r_ind] += jac * S.H[testx, testy, kx, ky] * S.H[trialx, trialy, kx, ky] * wgx[kx] * wgy[ky]
                end
            end
        end
    end

    #display(S.H)
    for i in bcs
        M[i, i] = 1.0
        A[i, i] = 1.0
    end


    #println(M)


    evals, efuncs = eigen(M, A)
    #evals, efuncs = eigen(Hermitian(M), Hermitian(A))

    #return evals, efuncs

    bc_evs = length(unique(bcs))

    #println(length(evals) - bc_evs)

    ef = zeros(length(evals)-bc_evs, Nx, Ny)

    ef = zeros(ComplexF64, length(evals), Nx, Ny)



    for i in 1:length(evals) #- bc_evs

        for j in 1:Nx * Ny

            #think perhaps this is not working...
            #seems to be, but resulting efuncs are kind of cooked
            #while evals are perf.
            xind, yind = index_to_grid_linear(j, Nx, Ny)
            ef[i, xind, yind] = efuncs[j, i]
        end
    end

        #for j in 1:Nx, k in 1:Ny

        #    ef[i, j, k] = efuncs[grid_to_index_linear(j, k, 1, 1, Nx, Ny), i+bc_evs]
        #end

    #println(M)

    display(bc_evs)

    return evals[1+bc_evs:end], ef[1+bc_evs:end, :, :]

end

function cubic(Nx, Ny, gpx, gpy)


    #xgrid = LinRange(0, π, Nx)#+1)[1:end-1]
    #ygrid = LinRange(0, π, Ny)#+1)[1:end-1]
    xgrid = LinRange(0, 1, Nx)#+1)[1:end-1]
    ygrid = LinRange(0, 2π, Ny+1)[1:end-1]

    ξx, wgx = gausslegendre(gpx)
    ξy, wgy = gausslegendre(gpy)

    S = hermite_basis(ξx, ξy)

    mat_size = cubic_mat_size(Nx, Ny)

    M = zeros(Float64, mat_size, mat_size)

    A = zeros(Float64, mat_size, mat_size)

    bcs = compute_cubic_bcs(Nx, Ny)
    

    for i in 1:Nx-1, j in 1:Ny

        x, y, dx, dy = linear_local_to_global(i, j, ξx, ξy, xgrid, ygrid)

        jac = dx * dy / 4

        for trialx in 1:4, trialy in 1:4

            r_ind = grid_to_index_cubic(i, j, trialx, trialy, Nx, Ny)

            if r_ind in bcs
                continue
            end

            for testx in 1:4, testy in 1:4

                l_ind = grid_to_index_cubic(i, j, testx, testy, Nx, Ny)

                if l_ind in bcs
                    continue
                end

                for kx in 1:gpx, ky in 1:gpy

                    #simple cartesian case
                    #M[l_ind, r_ind] += jac * (S.dHx1[testx, testy, kx, ky] * S.dHx1[trialx, trialy, kx, ky] * 2/dx * 2/dx + S.dHx2[testx, testy, kx, ky] * S.dHx2[trialx, trialy, kx, ky] * 2/dy * 2/dy) * wgx[kx] * wgy[ky]
                    #A[l_ind, r_ind] += jac * S.H[testx, testy, kx, ky] * S.H[trialx, trialy, kx, ky] * wgx[kx] * wgy[ky]

                    #cartesian case with artificail increase in derivatives.
                    #M[l_ind, r_ind] += jac * (S.ddHx1x2[testx, testy, kx, ky] * S.ddHx1x2[trialx, trialy, kx, ky] * (2/dx)^2 * (2/dy)^2 + S.ddHx2x2[testx, testy, kx, ky] * S.ddHx2x2[trialx, trialy, kx, ky] * (2/dy)^2 * (2/dy)^2) * wgx[kx] * wgy[ky]
                    #A[l_ind, r_ind] += jac * S.dHx2[testx, testy, kx, ky] * S.dHx2[trialx, trialy, kx, ky] * wgx[kx] * wgy[ky] * (2/dy)^2

                    #Polar, closer to ideal solution we want.
                    #seems to be working now, the degeneracy in n is kind of odd.
                    #M[l_ind, r_ind] += jac * (+S.dHx1[testx, testy, kx, ky] * S.dHx1[trialx, trialy, kx, ky] * (2/dx)^2 
                    #                          - S.H[testx, testy, kx, ky] / x[kx] * S.dHx1[trialx, trialy, kx, ky] * (2/dx)
                    #                          + S.dHx2[testx, testy, kx, ky] * S.dHx2[trialx, trialy, kx, ky] / x[kx]^2 * (2/dy)^2) * wgx[kx] * wgy[ky]

                    #A[l_ind, r_ind] += jac * S.H[testx, testy, kx, ky] * S.H[trialx, trialy, kx, ky] * wgx[kx] * wgy[ky]

                    #polar with extra double derivative of x2. aka θ
                    #this gives perfect eigenvalues
                    #but eigenfunctions are cooked af.
                    #mayhap this will work in 3d with z instead of θ?
                    M[l_ind, r_ind] += jac * (+S.ddHx1x2[testx, testy, kx, ky] * S.ddHx1x2[trialx, trialy, kx, ky] * (2/dx)^2 
                                              - S.dHx2[testx, testy, kx, ky] / x[kx] * S.ddHx1x2[trialx, trialy, kx, ky] * (2/dx)
                                              + S.ddHx2x2[testx, testy, kx, ky] * S.ddHx2x2[trialx, trialy, kx, ky] / x[kx]^2 * (2/dy)^2) * wgx[kx] * wgy[ky] * (2/dy)^2

                    A[l_ind, r_ind] += jac * S.dHx2[testx, testy, kx, ky] * S.dHx2[trialx, trialy, kx, ky] * wgx[kx] * wgy[ky] * (2/dy)^2

                end
            end
        end
    end

    #display(S.H)
    for i in bcs
        M[i, i] = 1.0
        A[i, i] = 1.0
    end


    #display(M)
    #println(M)


    #evals, efuncs = eigen(Hermitian(M), Hermitian(A))
    #evals, efuncs = eigen(M, A)
    nev = 10
    evals, efuncs = eigs(M, A, nev=nev, ritzvec=true, sigma=20)

    #return evals, efuncs

    #bc_evs = length(unique(bcs))

    #println(length(evals) - bc_evs)

    #ef = zeros(length(evals)-bc_evs, Nx, Ny)

    #ef = zeros(length(evals), Nx, Ny)
    #ef = zeros(ComplexF64, length(evals), Nx, Ny)
    ef = zeros(ComplexF64, nev, Nx, Ny)



    #for i in 1:length(evals) #- bc_evs
    for i in 1:nev

        for j in 1:4:4*Nx * Ny

            xind, yind = index_to_grid_cubic(j, Nx, Ny)
            #display((xind, yind))
            ef[i, xind, yind] = efuncs[j, i]
        end
    end


    #println(M)

    #println(bcs)

    return evals, ef
    #return evals[1+bc_evs:end], ef[1+bc_evs:end, :, :]

end
