
function cubic_mat_size(Nx::Int64, Ny::Int64, Nz::Int64)
    return 8*Nx*Ny*Nz
end

function compute_linear_bcs(Nx::Int64, Ny::Int64, Nz::Int64)
    leftx = 1:Ny
    rightx = (Nx-1)*Ny:Nx*Ny

    lefty = 1:Ny:Nx*Ny

    righty = Nx:Ny:Nx*Ny

    #[x1y1, x1y2, ... x1, yN, x2y1, x2y2, .... xNy1, xNy2, ..., xNyN]
    #return collect(rightx)

    return vcat(leftx, rightx)
    return vcat(leftx, lefty, rightx, righty)
end

function cubic_bcs_cart(Nx::Int64, Ny::Int64, Nz::Int64)


    bcs = Int64[]
    for i in 1:8*Nx*Ny*Nz
        xind, yind, zind, h = index_to_grid_cubic(i, Nx, Ny, Nz)

        if h==1 || h==2 || h==3 || h==4
            if xind == 1 || xind == Nx
                push!(bcs, i)
            end
        end
        if h==1 || h==2 || h==5 || h==6
            if yind == 1 || yind == Ny
                push!(bcs, i)
            end
        end
        if h==1 || h==3 || h==5 || h==7
            if zind == 1 || zind == Nz
                push!(bcs, i)
            end
        end
    end
    return bcs
end

    

function compute_cubic_bcs(Nx::Int64, Ny::Int64, Nz::Int64)
    Nx1 = Nx
    Nx2 = Ny
    Nx3 = Nz
    #note this only gets r boundaries, x2, x3 is assumed periodic and handled elsewhere.
    
    #8 basis functions are produced from the tensor product of 2 basis functions in 3d.
    #Boundary condition is applied to the first basis function in r
    #this occurs in the first 4/8 3d basis functions.
    left_boundary1 = 1:8:8*Nx2 * Nx3
    left_boundary2 = 2:8:8*Nx2 * Nx3
    left_boundary3 = 3:8:8*Nx2 * Nx3
    left_boundary4 = 4:8:8*Nx2 * Nx3
   
    right_boundary1 = 1 + (Nx1 - 1) * 8 * Nx2 * Nx3:8:8*Nx1 * Nx2 * Nx3
    right_boundary2 = 2 + (Nx1 - 1) * 8 * Nx2 * Nx3:8:8*Nx1 * Nx2 * Nx3
    right_boundary3 = 3 + (Nx1 - 1) * 8 * Nx2 * Nx3:8:8*Nx1 * Nx2 * Nx3
    right_boundary4 = 4 + (Nx1 - 1) * 8 * Nx2 * Nx3:8:8*Nx1 * Nx2 * Nx3

    return vcat(left_boundary1, left_boundary2, left_boundary3, left_boundary4, right_boundary1, right_boundary2, right_boundary3, right_boundary4)
end

function cubic(Nx, Ny, Nz, gpx, gpy, gpz)


    xgrid = LinRange(0, π, Nx)#+1)[1:end-1]
    ygrid = LinRange(0, π, Ny)#+1)[1:end-1]
    zgrid = LinRange(0, π, Nz)#+1)[1:end-1]
    xgrid = LinRange(0, 1, Nx)#+1)[1:end-1]
    ygrid = LinRange(0, 2π, Ny+1)[1:end-1]
    zgrid = LinRange(0, 2π, Nz+1)[1:end-1]

    ξx, wgx = gausslegendre(gpx)
    ξy, wgy = gausslegendre(gpy)
    ξz, wgz = gausslegendre(gpz)

    S = hermite_basis(ξx, ξy, ξz)

    mat_size = cubic_mat_size(Nx, Ny, Nz)

    M = zeros(Float64, mat_size, mat_size)

    A = zeros(Float64, mat_size, mat_size)

    bcs = compute_cubic_bcs(Nx, Ny, Nz)
    #bcs = cubic_bcs_cart(Nx, Ny, Nz)
    

    for i in 1:Nx-1, j in 1:Ny, k in 1:Nz
    #for i in 1:Nx-1, j in 1:Ny-1, k in 1:Nz-1

        x, y, z, dx, dy, dz = linear_local_to_global(i, j, k, ξx, ξy, ξz, xgrid, ygrid, zgrid)

        jac = dx * dy * dz / 8

        for trialx in 1:4, trialy in 1:4, trialz in 1:4

            r_ind = grid_to_index_cubic(i, j, k, trialx, trialy, trialz, Nx, Ny, Nz)

            if r_ind in bcs
                continue
            end

            for testx in 1:4, testy in 1:4, testz in 1:4

                l_ind = grid_to_index_cubic(i, j, k, testx, testy, testz, Nx, Ny, Nz)

                if l_ind in bcs
                    continue
                end

                for kx in 1:gpx, ky in 1:gpy, kz in 1:gpz

                    #simple cartesian case
                    #M[l_ind, r_ind] += jac * (S.dHx1[testx, testy, testz, kx, ky, kz] * S.dHx1[trialx, trialy, trialz, kx, ky, kz] * 2/dx * 2/dx 
                    #                          + S.dHx2[testx, testy, testz, kx, ky, kz] * S.dHx2[trialx, trialy, trialz, kx, ky, kz] * 2/dy * 2/dy
                    #                          + S.dHx3[testx, testy, testz, kx, ky, kz] * S.dHx3[trialx, trialy, trialz, kx, ky, kz] * (2/dz)^2) * wgx[kx] * wgy[ky] * wgz[kz]
                    #A[l_ind, r_ind] += jac * S.H[testx, testy, testz, kx, ky, kz] * S.H[trialx, trialy, trialz, kx, ky, kz] * wgx[kx] * wgy[ky] * wgz[kz]

                    #cartesian case with artificail increase in derivatives.
                    #touch less accurate, but this does still work
                    #looks like the accuracy is most noticible for case
                    #where kz is higher than others, as z is the most unstable.
                    #M[l_ind, r_ind] += jac * (S.ddHx1x3[testx, testy, testz, kx, ky, kz] * S.ddHx1x3[trialx, trialy, trialz, kx, ky, kz] * (2/dx)^2 * (2/dz)^2 
                    #                          + S.ddHx2x3[testx, testy, testz, kx, ky, kz] * S.ddHx2x3[trialx, trialy, trialz, kx, ky, kz] * (2/dy)^2 * (2/dz)^2
                    #                          + S.ddHx3x3[testx, testy, testz, kx, ky, kz] * S.ddHx3x3[trialx, trialy, trialz, kx, ky, kz] * (2/dz)^4) * wgx[kx] * wgy[ky] * wgz[kz]
                    #A[l_ind, r_ind] += jac * S.dHx3[testx, testy, testz, kx, ky, kz] * S.dHx3[trialx, trialy, trialz, kx, ky, kz] * (2/dz)^2 * wgx[kx] * wgy[ky] * wgz[kz]

                    #Polar, closer to ideal solution we want.
                    #seems to be working now, the degeneracy in n is kind of odd.
                    #M[l_ind, r_ind] += jac * (+S.dHx1[testx, testy, testz, kx, ky, kz] * S.dHx1[trialx, trialy, trialz, kx, ky, kz] * (2/dx)^2 
                    #                          - S.H[testx, testy, testz, kx, ky, kz] / x[kx] * S.dHx1[trialx, trialy, trialz, kx, ky, kz] * (2/dx)
                    #                          + S.dHx2[testx, testy, testz, kx, ky, kz] * S.dHx2[trialx, trialy, trialz, kx, ky, kz] / x[kx]^2 * (2/dy)^2
                    #                          + S.dHx3[testx, testy, testz, kx, ky, kz] * S.dHx3[trialx, trialy, trialz, kx, ky, kz] * (2/dz)^2) * wgx[kx] * wgy[ky] * wgz[kz]

                    #A[l_ind, r_ind] += jac * S.H[testx, testy, testz, kx, ky, kz] * S.H[trialx, trialy, trialz, kx, ky, kz] * wgx[kx] * wgy[ky] * wgz[kz]

                    #polar with extra double derivative of x3. aka ζ
                    M[l_ind, r_ind] += jac * (+S.ddHx1x3[testx, testy, testz, kx, ky, kz] * S.ddHx1x3[trialx, trialy, trialz, kx, ky, kz] * (2/dx)^2 *(2/dz)^2
                                              - S.dHx3[testx, testy, testz, kx, ky, kz] / x[kx] * S.ddHx1x3[trialx, trialy, trialz, kx, ky, kz] * (2/dx) * (2/dz)^2
                                              + S.ddHx2x3[testx, testy, testz, kx, ky, kz] * S.ddHx2x3[trialx, trialy, trialz, kx, ky, kz] / x[kx]^2 * (2/dy)^2 * (2/dz)^2
                                              + S.ddHx3x3[testx, testy, testz, kx, ky, kz] * S.ddHx3x3[trialx, trialy, trialz, kx, ky, kz] * (2/dz)^4) * wgx[kx] * wgy[ky] * wgz[kz]

                    A[l_ind, r_ind] += jac * S.dHx3[testx, testy, testz, kx, ky, kz] * S.dHx3[trialx, trialy, trialz, kx, ky, kz] * wgx[kx] * wgy[ky] * wgz[kz] * (2/dz)^2

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
    #evals, efuncs = eigs(M, A, nev=nev, ritzvec=true, sigma=7)

    evals, efuncs = eigen(M, A)
    nev = length(evals)
    #return evals, efuncs

    bc_evs = length(unique(bcs))

    #println(length(evals) - bc_evs)

    #ef = zeros(length(evals)-bc_evs, Nx, Ny)

    #ef = zeros(length(evals), Nx, Ny)
    #ef = zeros(ComplexF64, length(evals), Nx, Ny)
    ef = zeros(ComplexF64, nev, Nx, Ny, Nz)



    #for i in 1:length(evals) #- bc_evs
    for i in 1:nev

        for j in 1:8:8*Nx * Ny*Nz

            xind, yind, zind, _ = index_to_grid_cubic(j, Nx, Ny, Nz)
            #display((xind, yind))
            ef[i, xind, yind, zind] = efuncs[j, i]
        end
    end


    #println(M)

    display(bc_evs)

    #return evals, ef
    return evals[1+bc_evs:end], ef[1+bc_evs:end, :, :, :]

end
