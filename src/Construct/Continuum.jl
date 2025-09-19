""" 
    continuum(prob::ProblemT, grids::ContGridsT)

Finds the continuum by solving the second order derivatives on each flux surface individually. Much faster than reconstructing the continuum from the full solver, but cannot handle islands, resistivity and won't find any global modes. This function assumes the spectral method is used in x2 and x3.
"""
function continuum(prob::ProblemT, grids::ContGridsT)


    x1grid, x2grid, x3grid = inst_grids(grids)

    Nx2 = length(x2grid)
    Nx3 = length(x3grid)

    mlist = mode_list(grids.x2)
    nlist = mode_list(grids.x3)


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
    I = local_matrix_size(grids)
    W = local_matrix_size(grids)

    Imat = zeros(ComplexF64, mat_size, mat_size)
    Wmat = zeros(ComplexF64, mat_size, mat_size)

    #sub matrices storing the relavant parts for continuum computation.
    Icont = zeros(1, 1, 1, Nx2, Nx3)
    Wcont = zeros(3, 3, 1, Nx2, Nx3)

    #Plans for efficient fourier transform
    pI = plan_fft(Icont, [4, 5])
    pW = plan_fft(Wcont, [4, 5])


    #struct for storing temp matrices used for the weakform.
    tm = TM()

    #now we loop through the grid
    for (i, x1) in enumerate(x1grid) 

        #computes the full weakform
        W_and_I!(W, I, B, met, prob, [x1], x2grid, x3grid, tm)

        #For the continuum we extract the relevant second derivative parts
        Icont = pI * I[[1], [1], :, :, :]
        Wcont = pW * W[[1, 5, 6], [1, 5, 6], :, :, :]
        
        Imat .= 0.0 + 0.0im
        Wmat .= 0.0 + 0.0im

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
                Imat[left_ind, right_ind] += Icont[1, 1, 1, mind, nind]

                #contract over the different derivatives.
                for j in 1:3
                    for k in 1:3
                        Wmat[left_ind, right_ind] += Ψ[j] * Wcont[j, k, 1, mind, nind] * Φ[k]
                    end
                         
                end
            end
        end

        
        vals = eigvals(Hermitian(Wmat), Hermitian(Imat))

        #normalise the eigenvalues
        #abs here is not ideal.
        ωlist[i, :] = prob.geo.R0 * sqrt.(abs.(vals))

    end

    #don't think modelabs will actually work for this case but whatever.
    #actually, surely it does?
    #nah they are returned sorted.

    evals_ω = ComplexF64[]
    evals_x1 = Float64[]
    evals_ml = Tuple{Int64, Int64}[]

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
    I = local_matrix_size(grids)
    W = local_matrix_size(grids)

    Imat = zeros(ComplexF64, mat_size, mat_size)
    Wmat = zeros(ComplexF64, mat_size, mat_size)

    #so turns out our original continuum funciton is garbage.
    Icont = zeros(1, 1, 1, Nx2, Nx3)

    Wcont = zeros(3, 3, 1, Nx2, Nx3)

    #probably not how this works.
    pI = plan_fft(Icont, [4, 5])
    pW = plan_fft(Wcont, [4, 5])

    #structs storing temporary arrays.
    CT = CoordTransformT()
    tm = TM()

    for (i, s) in enumerate(x1grid)


        #computes the weakform and coordinate transform
        W_and_I!(W, I, tor_B, tor_met, qfm_B, qfm_met, prob, [s], x2grid, x3grid, tm, surf_itp, CT, sd)


        Icont = pI * I[[1], [1], :, :, :]
        Wcont = pW * W[[1, 5, 6], [1, 5, 6], :, :, :]
        

        Imat .= 0.0
        Wmat .= 0.0

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
                Imat[left_ind, right_ind] += Icont[1, 1, 1, mind, nind]

                #contract over the different derivatives.
                for j in 1:3
                    for k in 1:3
                        Wmat[left_ind, right_ind] += Ψ[j] * Wcont[j, k, 1, mind, nind] * Φ[k]
                    end
                         
                end
            end
        end

        #not sure if this will always be Hermitian
        vals = eigvals(Hermitian(Wmat), Hermitian(Imat))
        #vals = eigvals(Wmat, Imat)

        ωlist[i, :] = prob.geo.R0 * sqrt.(abs.(vals))
    end

    #should probably return an evals Obj.
    return ωlist

end
