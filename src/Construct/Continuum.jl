""" 
    continuum(prob::ProblemT, grids::ContGridsT)

Finds the continuum by solving the second order derivatives on each flux surface individually.
Much faster than reconstructing the continuum from the full solver.
This function assumes the spectral method is used in x2 and x3.
"""
function continuum(prob::ProblemT, grids::ContGridsT)

    #instantiate the grids into arrays.
    x1grid, x2grid, x3grid = inst_grids(grids)

    #for spectral method we need the length of the arrays
    Nx2 = length(x2grid)
    Nx3 = length(x3grid)

    #and the list of modes to consider.
    mlist = mode_list(grids.x2)
    nlist = mode_list(grids.x3)

    #initialise the two structs to store the metric and the magnetic field.
    met = MetT()
    B = BFieldT()

    #size of matrices
    mat_size = matrix_size(grids)

    #array to store the eigenvalues
    ωlist = zeros(ComplexF64, grids.x1.N, mat_size)

    #generalised eval problem PΦ = ω^2 Q Φ
    #these matrices store the local contribution, i.e. at each grid point, for the global matrices Q and P.
    Q = init_local_matrix(grids)
    P = init_local_matrix(grids)

    #global matrices.
    Qmat = zeros(ComplexF64, mat_size, mat_size)
    Pmat = zeros(ComplexF64, mat_size, mat_size)

    #sub matrices storing the relavant parts for continuum computation.
    Qcont = zeros(ComplexF64, 1, 1, 1, Nx2, Nx3)
    Pcont = zeros(ComplexF64, 3, 3, 1, Nx2, Nx3)

    #Plans for efficient fourier transform
    pQ = plan_fft(Qcont, [4, 5])
    pP = plan_fft(Pcont, [4, 5])

    #struct for storing temp matrices used for the weakform.
    tm = TM()

    #now we loop through the grid
    for (i, x1) in enumerate(x1grid) 

        #computes the full weakform
        weak_form!(P, Q, B, met, prob, [x1], x2grid, x3grid, tm)

        #For the continuum we extract the relevant second derivative parts
        Qcont .= pQ * Q[[1], [1], :, :, :]
        Pcont .= pP * P[[1, 5, 6], [1, 5, 6], :, :, :]
        
        Qmat .= 0.0 + 0.0im
        Pmat .= 0.0 + 0.0im

        #loop over the fourier components of the trial function
        for (k1,m1) in enumerate(mlist), (l1, n1) in enumerate(nlist)

            left_ind = grid_to_index(k1, l1, grids)

            #equivalent to create local basis. Much simpler in continuum case.
            Φ = [1, m1*1im, n1*1im]
            
            for (k2, m2) in enumerate(mlist), (l2, n2) in enumerate(nlist)

                right_ind = grid_to_index(k2, l2, grids)
                
                Ψ = [1, -m2*1im, -n2*1im]

                mind = mod(k1-k2 + Nx2, Nx2) + 1
                nind = mod(l1-l2 + Nx3, Nx3) + 1

                #silly matrix dimensions, but keeps it consistent.
                Qmat[left_ind, right_ind] += Qcont[1, 1, 1, mind, nind]

                #contract over the different derivatives.
                for j in 1:3
                    for k in 1:3
                        Pmat[left_ind, right_ind] += Ψ[j] * Pcont[j, k, 1, mind, nind] * Φ[k]
                    end
                end
            end
        end

        #eigenvalues are found for each radial point.
        vals = eigvals(Hermitian(Pmat), Hermitian(Qmat))

        ωlist[i, :] = vals

    end

    return ωlist

end



"""
    qfm_continuum(prob::ProblemT, grids::ContGridsT, surfs::Array{QFMSurfaceT})

Finds the continuum by solving the second order derivatives on each flux surface individually using qfm surface.
This function assumes the spectral method is used in x2 and x3.
"""
function continuum(prob::ProblemT, grids::ContGridsT, surfs::Array{QFMSurfaceT})

    #instantiate the grids into arrays.
    x1grid, x2grid, x3grid = inst_grids(grids)

    #for spectral method we need the length of the arrays
    Nx2 = length(x2grid)
    Nx3 = length(x3grid)

    #and the list of modes to consider.
    mlist = mode_list(grids.x2)
    nlist = mode_list(grids.x3)

    #creates the interpolation from the surfaces
    surf_itp, sd = create_surf_itp(surfs);

    #define structures to store the old and new metric and magnetic field
    tor_met = MetT()
    qfm_met = MetT()
    tor_B = BFieldT()
    qfm_B = BFieldT()

    #size of matrices
    mat_size = matrix_size(grids)

    #array to store the eigenvalues
    ωlist = zeros(ComplexF64, grids.x1.N, mat_size)

    #generalised eval problem PΦ = ω^2 Q Φ
    #these matrices store the local contribution, i.e. at each grid point, for the global matrices Q and P.
    Q = init_local_matrix(grids)
    P = init_local_matrix(grids)

    #global matrices.
    Qmat = zeros(ComplexF64, mat_size, mat_size)
    Pmat = zeros(ComplexF64, mat_size, mat_size)

    #sub matrices storing the relavant parts for continuum computation.
    Qcont = zeros(ComplexF64, 1, 1, 1, Nx2, Nx3)
    Pcont = zeros(ComplexF64, 3, 3, 1, Nx2, Nx3)

    #Plans for efficient fourier transform
    pQ = plan_fft(Qcont, [4, 5])
    pP = plan_fft(Pcont, [4, 5])

    #struct for storing temp matrices used for the weakform.
    tm = TM()

    #struct for storing the intermediate data for the coordinate transform
    CT = CoordTransformT()

    #now we loop through the grid
    for (i, x1) in enumerate(x1grid) 

        #computes the weakform and coordinate transform
        weak_form!(P, Q, tor_B, tor_met, qfm_B, qfm_met, prob, [x1], x2grid, x3grid, tm, surf_itp, CT, sd)

        #For the continuum we extract the relevant second derivative parts
        Qcont .= pQ * Q[[1], [1], :, :, :]
        Pcont .= pP * P[[1, 5, 6], [1, 5, 6], :, :, :]
        
        Qmat .= 0.0 + 0.0im
        Pmat .= 0.0 + 0.0im

        #loop over the fourier components of the trial function
        for (k1,m1) in enumerate(mlist), (l1, n1) in enumerate(nlist)

            left_ind = grid_to_index(k1, l1, grids)

            #equivalent to create local basis. Much simpler in continuum case.
            Φ = [1, m1*1im, n1*1im]
            
            for (k2, m2) in enumerate(mlist), (l2, n2) in enumerate(nlist)

                right_ind = grid_to_index(k2, l2, grids)
                
                Ψ = [1, -m2*1im, -n2*1im]

                mind = mod(k1-k2 + Nx2, Nx2) + 1
                nind = mod(l1-l2 + Nx3, Nx3) + 1

                #silly matrix dimensions, but keeps it consistent.
                Qmat[left_ind, right_ind] += Qcont[1, 1, 1, mind, nind]

                #contract over the different derivatives.
                for j in 1:3
                    for k in 1:3
                        Pmat[left_ind, right_ind] += Ψ[j] * Pcont[j, k, 1, mind, nind] * Φ[k]
                    end
                end
            end
        end

        #eigenvalues are found for each radial point.
        vals = eigvals(Hermitian(Pmat), Hermitian(Qmat))

        ωlist[i, :] = vals

    end

    return ωlist

end

