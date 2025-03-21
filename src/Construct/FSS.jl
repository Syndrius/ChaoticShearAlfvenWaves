"""
    construct(prob::ProblemT, grids::FSSGridsT)

Constructs the two matrices using the WeakForm of the SAW governing equation. Uses Finite elements with cubic Hermite polynomials in r, and the fourier spectral method in θ and ζ. Returns two sparse matrices.

### Args
prob::ProblemT - Struct containing the functions and parameters that define the problem we are solving
grids::FSSGridT - Grids to solve over.
"""
function construct(prob::ProblemT, grids::FSSGridsT)

    #instantiate the grids into arrays.
    rgrid, θgrid, ζgrid = inst_grids(grids)

    #for spectral method we need the length of the arrays
    Nθ = length(θgrid)
    Nζ = length(ζgrid)

    #and the list of modes to consider.
    mlist = mode_list(grids.θ)
    nlist = mode_list(grids.ζ)

    #initialise the two structs to store the metric and the magnetic field.
    met = MetT()
    B = BFieldT()

    #compute the gaussian qudrature points for finite elements.
    ξ, wg = gausslegendre(grids.r.gp) 

    #Gets the Hermite basis for the radial grid.
    S = hermite_basis(ξ)

    
    #creates the trial and test function arrays.
    #these store the basis functions for each derivative
    #and finite elements basis 
    Φ = init_basis_function(grids)
    Ψ = init_basis_function(grids)
    
    #arrays to store the row, column and data of each matrix element
    #used for constructing sparse matrices.
    rows = Array{Int64}(undef, 0) 
    cols = Array{Int64}(undef, 0) 
    Idata = Array{ComplexF64}(undef, 0)
    Wdata = Array{ComplexF64}(undef, 0)


    #determines the indicies in the matrices corresponding to the dirichlet boundary conditions for r.
    boundary_inds = compute_boundary_inds(grids)
    

    #generalised eval problem WΦ = ω^2 I Φ
    #these matrices store the local contribution, i.e. at each grid point, for the global matrices I and W.
    I = local_matrix_size(grids)
    W = local_matrix_size(grids)


    #creates a fft plan for efficient fft used in spectral method.
    p = plan_fft!(W, [4, 5])

    #initialises a struct storing temporary matrices used in the weak form.
    tm = TM()
   
    #main loop
    for i in 1:grids.r.N-1

        #takes the local ξ array to a global r array around the grid point.
        r, dr = local_to_global(i, ξ, rgrid)

        #jacobian of the local to global transformation.
        jac = dr/2 

        #computes the contribution to the W and I matrices.
        W_and_I!(W, I, met, B, prob, r, θgrid, ζgrid, tm)
        
        #fft the two matrices.
        p * W
        p * I

        #loop over the fourier components of the trial function
        for (k1,m1) in enumerate(mlist), (l1, n1) in enumerate(nlist)

            #adjust the basis functions to the current coordinates/mode numbers considered.
            create_local_basis!(Φ, S, m1, n1, jac)

            for (k2, m2) in enumerate(mlist), (l2, n2) in enumerate(nlist)

                #negatives for conjugate in the test function.
                create_local_basis!(Ψ, S, -m2, -n2, jac)

                #extract the relevant indicies from the fft'ed matrices.
                mind = mod(k1-k2 + Nθ, Nθ) + 1
                nind = mod(l1-l2 + Nζ, Nζ) + 1

                #loop over the Hermite elements for the trial function
                for trialsf in 1:4

                    #determines the matrix index for the trial function
                    right_ind = grid_to_index(i, k1, l1, trialsf, grids)

                    #and for the test function.
                    for testsf in 1:4

                        #index for the test function
                        left_ind = grid_to_index(i, k2, l2, testsf, grids)

                        #only check for boundaries if this is true
                        #no other i's can possibly give boundaries
                        if i==1 || i==grids.r.N-1

                            if left_ind == right_ind && left_ind in boundary_inds

                                #diagonals for boundary conditions are set to 1.
                                push!(rows, left_ind)
                                push!(cols, right_ind)
                                push!(Wdata, 1.0 + 0.0im)
                                push!(Idata, 1.0 + 0.0im)
                                                 
                            #otherwise the boundaries are set to zero, which for sparse matrices
                            #is the same as leaving blank.
                            elseif left_ind in boundary_inds
                                continue
                            elseif right_ind in boundary_inds
                                continue

                            #otherwise a regular case for these indicies.
                            else
                                
                                #integrate the local contribution to our matrices.
                                Wsum = @views gauss_integrate(Ψ[testsf, :, :], Φ[trialsf, :, :], W[:, :, :, mind, nind], wg, jac, grids.r.gp)
                                Isum = @views gauss_integrate(Ψ[testsf, :, :], Φ[trialsf, :, :], I[:, :, :, mind, nind], wg, jac, grids.r.gp)

                                #adds the local contribution to the global structure
                                push!(Wdata, Wsum)
                                push!(Idata, Isum)
                                push!(rows, left_ind)
                                push!(cols, right_ind)
                            end
                        else

                            #integrate the local contribution to our matrices.
                            Wsum = @views gauss_integrate(Ψ[testsf, :, :], Φ[trialsf, :, :], W[:, :, :, mind, nind], wg, jac, grids.r.gp)
                            Isum = @views gauss_integrate(Ψ[testsf, :, :], Φ[trialsf, :, :], I[:, :, :, mind, nind], wg, jac, grids.r.gp)

                            #adds the local contribution to the global structure
                            push!(Wdata, Wsum)
                            push!(Idata, Isum)
                            push!(rows, left_ind)
                            push!(cols, right_ind) 
                        
                        end
                    end
                end

            end

        end
    end

    #construct the sparse matrix.
    Wmat = sparse(rows, cols, Wdata)
    Imat = sparse(rows, cols, Idata)
    
    #display(Matrix(Imat))

    return Wmat, Imat
end




function construct(prob::ProblemT, grids::FSSGridsT, surfs::Array{QFMSurfaceT})

    #instantiate the grids into arrays.
    #note that the inputs are the new coords.
    sgrid, ϑgrid, φgrid = inst_grids(grids)

    #for spectral method we need the length of the arrays
    Nϑ = length(ϑgrid)
    Nφ = length(φgrid)

    #and the list of modes to consider.
    mlist = mode_list(grids.θ)
    nlist = mode_list(grids.ζ)

    #initialise the two structs to store the metric and the magnetic field.
    #one for each coordinate system
    tor_met = MetT()
    qfm_met = MetT()
    tor_B = BFieldT()
    qfm_B = BFieldT()

    #creates the interpolations for the surfaces.
    surf_itp = create_surf_itp(surfs)

    #compute the gaussian qudrature points for finite elements.
    ξ, wg = gausslegendre(grids.r.gp) 

    #Gets the Hermite basis for the radial grid.
    S = hermite_basis(ξ)

    
    #creates the trial and test function arrays.
    #these store the basis functions for each derivative
    #and finite elements basis 
    Φ = init_basis_function(grids)
    Ψ = init_basis_function(grids)
    
    #arrays to store the row, column and data of each matrix element
    #used for constructing sparse matrices.
    rows = Array{Int64}(undef, 0) 
    cols = Array{Int64}(undef, 0) 
    Idata = Array{ComplexF64}(undef, 0)
    Wdata = Array{ComplexF64}(undef, 0)


    #determines the indicies in the matrices corresponding to the dirichlet boundary conditions for r.
    boundary_inds = compute_boundary_inds(grids)
    

    #generalised eval problem WΦ = ω^2 I Φ
    #these matrices store the local contribution, i.e. at each grid point, for the global matrices I and W.
    I = local_matrix_size(grids)
    W = local_matrix_size(grids)


    #creates a fft plan for efficient fft used in spectral method.
    p = plan_fft!(W, [4, 5])


    #initialises a struct storing temporary matrices used in the weak form.
    tm = TM()

    #struct for storing the intermediate data for the coordinate transform
    CT = CoordTsfmT()
   
    #main loop
    for i in 1:grids.r.N-1

        #takes the local ξ array to a global r array around the grid point.
        s, ds = local_to_global(i, ξ, sgrid)

        #jacobian of the local to global transformation.
        jac = ds/2 

        #computes the contribution to the W and I matrices.
        W_and_I!(W, I, tor_met, tor_B, qfm_met, qfm_B, prob, s, ϑgrid, φgrid, tm, surf_itp, CT)

        #display(W)
        
        #fft the two matrices.
        p * W
        p * I

        #loop over the fourier components of the trial function
        for (k1,m1) in enumerate(mlist), (l1, n1) in enumerate(nlist)

            #adjust the basis functions to the current coordinates/mode numbers considered.
            create_local_basis!(Φ, S, m1, n1, jac)

            for (k2, m2) in enumerate(mlist), (l2, n2) in enumerate(nlist)

                #negatives for conjugate in the test function.
                create_local_basis!(Ψ, S, -m2, -n2, jac)

                #extract the relevant indicies from the fft'ed matrices.
                mind = mod(k1-k2 + Nϑ, Nϑ) + 1
                nind = mod(l1-l2 + Nφ, Nφ) + 1

                #loop over the Hermite elements for the trial function
                for trialsf in 1:4

                    #determines the matrix index for the trial function
                    right_ind = grid_to_index(i, k1, l1, trialsf, grids)

                    #and for the test function.
                    for testsf in 1:4

                        #index for the test function
                        left_ind = grid_to_index(i, k2, l2, testsf, grids)

                        #only check for boundaries if this is true
                        #no other i's can possibly give boundaries
                        if i==1 || i==grids.r.N-1

                            if left_ind == right_ind && left_ind in boundary_inds

                                #diagonals for boundary conditions are set to 1.
                                push!(rows, left_ind)
                                push!(cols, right_ind)
                                push!(Wdata, 1.0 + 0.0im)
                                push!(Idata, 1.0 + 0.0im)
                                                 
                            #otherwise the boundaries are set to zero, which for sparse matrices
                            #is the same as leaving blank.
                            elseif left_ind in boundary_inds
                                continue
                            elseif right_ind in boundary_inds
                                continue

                            #otherwise a regular case for these indicies.
                            else
                                
                                #integrate the local contribution to our matrices.
                                Wsum = @views gauss_integrate(Ψ[testsf, :, :], Φ[trialsf, :, :], W[:, :, :, mind, nind], wg, jac, grids.r.gp)
                                Isum = @views gauss_integrate(Ψ[testsf, :, :], Φ[trialsf, :, :], I[:, :, :, mind, nind], wg, jac, grids.r.gp)

                                #adds the local contribution to the global structure
                                push!(Wdata, Wsum)
                                push!(Idata, Isum)
                                push!(rows, left_ind)
                                push!(cols, right_ind)
                            end
                        else

                            #integrate the local contribution to our matrices.
                            Wsum = @views gauss_integrate(Ψ[testsf, :, :], Φ[trialsf, :, :], W[:, :, :, mind, nind], wg, jac, grids.r.gp)
                            Isum = @views gauss_integrate(Ψ[testsf, :, :], Φ[trialsf, :, :], I[:, :, :, mind, nind], wg, jac, grids.r.gp)

                            #adds the local contribution to the global structure
                            push!(Wdata, Wsum)
                            push!(Idata, Isum)
                            push!(rows, left_ind)
                            push!(cols, right_ind) 
                        
                        end
                    end
                end

            end

        end
    end

    #construct the sparse matrix.
    Wmat = sparse(rows, cols, Wdata)
    Imat = sparse(rows, cols, Idata)

    #display(Matrix(Imat))

    return Wmat, Imat
end
