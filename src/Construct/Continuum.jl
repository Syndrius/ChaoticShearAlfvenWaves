""" 
    continuum(prob::ProblemT, grids::ContGridsT)

Finds the continuum by solving the second order derivatives on each flux surface individually. Much faster than reconstructing the continuum from the full solver, but cannot handle islands, resistivity and won't find any global modes. This function assumes the spectral method is used in x2 and x3.
"""
function continuum(prob::ProblemT, grids::ContGridsT)


    x1grid, x2grid, x3grid = inst_grids(grids)


    #seems to make no difference.
    #x3grid = range(0, 2π * prob.isls[1].m0, grids.x3.N * grids.x3.f_quad+1)[1:end-1]

    #display(x3grid)
    Nx2 = length(x2grid)
    Nx3 = length(x3grid)

    mlist = mode_list(grids.x2)
    nlist = mode_list(grids.x3)
    #nlist = grids.x3.start:grids.x3.incr*prob.isls[1].m0:grids.x3.stop*prob.isls[1].m0

    #display(collect(nlist))


    #this condition no longer works with the new setup.
    #probably need to check that there is only a single, A=0.0 island.
    #=
    if !isempty(prob.isls)
        display("Continuum can't handle islands you goose.")
        #exit()
    end
    =#

    if prob.flr.δ != 0.0
        display("Continuum can't handle resisitivity you goose.")
        #exit()
    end


    met = MetT()
    B = BFieldT()


    mat_size = matrix_size(grids)

    #array to store the eigenvalues
    ωlist = zeros(grids.x1.N, mat_size)

    #matrices for the local contribution to the global matrices.#
    Q = init_local_matrix(grids)
    P = init_local_matrix(grids)

    #should be a function!
    Qmat = zeros(ComplexF64, mat_size, mat_size)
    Pmat = zeros(ComplexF64, mat_size, mat_size)

    #sub matrices storing the relavant parts for continuum computation.
    Qcont = zeros(1, 1, 1, Nx2, Nx3)
    Pcont = zeros(3, 3, 1, Nx2, Nx3)

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
        #seems like this could be done in place tbh!
        Qcont = pQ * Q[[1], [1], :, :, :]
        Pcont = pP * P[[1, 5, 6], [1, 5, 6], :, :, :]
        
        Qmat .= 0.0 + 0.0im
        Pmat .= 0.0 + 0.0im

        #loop over the fourier components of the trial function
        for (k1,m1) in enumerate(mlist), (l1, n1) in enumerate(nlist)

            left_ind = grid_to_index(k1, l1, grids)

            #equivalent to create local basis. Much simpler in contninuum case.
            Φ = [1, m1*1im, n1*1im]
            
            for (k2, m2) in enumerate(mlist), (l2, n2) in enumerate(nlist)

                right_ind = grid_to_index(k2, l2, grids)
                
                Ψ = [1, -m2*1im, -n2*1im]

                mind = mod(k1-k2 + Nx2, Nx2) + 1
                nind = mod(l1-l2 + Nx3, Nx3) + 1

                #silly matrix dims, but keeps it consistent.
                Qmat[left_ind, right_ind] += Qcont[1, 1, 1, mind, nind]

                #contract over the different derivatives.
                for j in 1:3
                    for k in 1:3
                        Pmat[left_ind, right_ind] += Ψ[j] * Pcont[j, k, 1, mind, nind] * Φ[k]
                    end
                         
                end
            end
        end

        
        vals = eigvals(Hermitian(Pmat), Hermitian(Qmat))

        #normalise the eigenvalues
        #abs here is not ideal.
        #ωlist[i, :] = prob.geo.R0 * sqrt.(abs.(vals))
        #normalising is done elsewhere.
        ωlist[i, :] = vals

    end

    #don't think modelabs will actually work for this case but whatever.
    #actually, surely it does?
    #nah they are returned sorted.

    return ωlist

    #this needs to be moved to post-processing!
    #=
    evals_ω = ComplexF64[]
    evals_x1 = Float64[]
    evals_ml = Tuple{Qnt64, Qnt64}[]

    for i in 1:grids.x1.N
        for j in 1:grids.x2.N * grids.x3.N
            push!(evals_ω, ωlist[i, j])
            push!(evals_x1, x1grid[i])
            #don't know what these are!
            push!(evals_ml, (0, 0))
        end
    end

    #ideally this would return an evals like structure instead of this
    #TODO, either make this return evalsT, or create postprocessing step that does this
    #would help with making the plotting less shit, maybe trickier than thought though

    return EvalsT(evals_ω, evals_x1, evals_ml)
    =#
end




"""
    qfm_continuum(prob::ProblemT, grids::ContGridsT, surfs::Array{QFMSurfaceT})

Computes the continuum 
Finds the continuum by solving the second order derivatives on each flux surface individually using qfm surface. This function assumes the spectral method is used in x2 and x3.
"""
function continuum(prob::ProblemT, grids::ContGridsT, surfs::Array{QFMSurfaceT})
    
    #note that this should be in terms of (s, ϑ, φ) which is converted to (r, x2, x3) during the qfm process.

    x1grid, x2grid, x3grid = inst_grids(grids)

    Nx2 = length(x2grid)
    Nx3 = length(x3grid)

    mlist = mode_list(grids.x2)
    nlist = mode_list(grids.x3)

    #creates the interpolation from the surfaces
    surf_itp, sd = create_surf_itp(surfs);

    #define structures to store the old and new metric and magnetic field
    tor_met = MetT()
    qfm_met = MetT()
    tor_B = BFieldT()
    qfm_B = BFieldT()

    mat_size = matrix_size(grids)

    #array to store the results
    ωlist = zeros(grids.x1.N, mat_size)

    #pretty sure we have functions to do a bunch of this shite.
    Q = local_matrix_size(grids)
    P = local_matrix_size(grids)

    Qmat = zeros(ComplexF64, mat_size, mat_size)
    Pmat = zeros(ComplexF64, mat_size, mat_size)

    #so turns out our original continuum funciton is garbage.
    Qcont = zeros(1, 1, 1, Nx2, Nx3)

    Pcont = zeros(3, 3, 1, Nx2, Nx3)

    #probably not how this works.
    pQ = plan_fft(Qcont, [4, 5])
    pP = plan_fft(Pcont, [4, 5])

    #structs storing temporary arrays.
    CT = CoordTransformT()
    tm = TM()

    for (i, s) in enumerate(x1grid)


        #computes the weakform and coordinate transform
        weak_form!(P, Q, tor_B, tor_met, qfm_B, qfm_met, prob, [s], x2grid, x3grid, tm, surf_itp, CT, sd)


        Qcont = pQ * Q[[1], [1], :, :, :]
        Pcont = pP * P[[1, 5, 6], [1, 5, 6], :, :, :]
        

        Qmat .= 0.0
        Pmat .= 0.0

        for (k1,m1) in enumerate(mlist), (l1, n1) in enumerate(nlist)

            left_ind = grid_to_index(k1, l1, grids)

            #equivalent to create local basis. Much simpler in contninuum case.
            Φ = [1, m1*1im, n1*1im]
            
            for (k2, m2) in enumerate(mlist), (l2, n2) in enumerate(nlist)

                right_ind = grid_to_index(k2, l2, grids)
                
                Ψ = [1, -m2*1im, -n2*1im]

                mind = mod(k1-k2 + Nx2, Nx2) + 1
                nind = mod(l1-l2 + Nx3, Nx3) + 1

                #silly matrix dims, but keeps it consistent.
                Qmat[left_ind, right_ind] += Qcont[1, 1, 1, mind, nind]

                #contract over the different derivatives.
                for j in 1:3
                    for k in 1:3
                        Pmat[left_ind, right_ind] += Ψ[j] * Pcont[j, k, 1, mind, nind] * Φ[k]
                    end
                         
                end
            end
        end

        #not sure if this will always be Hermitian
        vals = eigvals(Hermitian(Pmat), Hermitian(Qmat))
        #vals = eigvals(Pmat, Qmat)

        ωlist[i, :] = prob.geo.R0 * sqrt.(abs.(vals))
    end

    #should probably return an evals Obj.
    return ωlist

end
