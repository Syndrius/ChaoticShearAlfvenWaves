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
    met = init_empty_met()
    B = init_empty_B()

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

    return Wmat, Imat
end



"""
    construct(prob::ProblemT, grids::FFSGridsT)

Constructs the two matrices using the WeakForm of the SAW governing equation. Uses Finite elements with cubic Hermite polynomials in r and θ and the fourier spectral method in ζ. Returns two sparse matrices.

### Args
prob::ProblemT - Struct containing the functions and parameters that define the problem we are solving
grids::FFSGridT - Grids to solve over.
"""
function construct(prob::ProblemT, grids::FFSGridsT)
    
    #instantiate the grids into arrays.
    rgrid, θgrid, ζgrid = inst_grids(grids)

    #for spectral method we need the length of the array
    Nζ = length(ζgrid)
    #and the mode list.
    nlist = mode_list(grids.ζ)


    #initialise the two structs to store the metric and the magnetic field.
    met = MetT()
    B = BFieldT()

    #compute the gaussian qudrature points for finite elements.
    ξr, wgr = gausslegendre(grids.r.gp) 
    ξθ, wgθ = gausslegendre(grids.θ.gp)

    
    #Gets the Hermite basis for the radial grid and poloidal grid.
    S = hermite_basis(ξr, ξθ)

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
    p = plan_fft!(W, [5])

    #initialises a struct storing temporary matrices used in the weak form.
    tm = TM()

    
    #main loop, θ goes to N for periodicity.
    for i in 1:grids.r.N-1, j in 1:grids.θ.N 

        #takes the local ξ arrays to a global arrays around the grid point.
        r, θ, dr, dθ = local_to_global(i, j, ξr, ξθ, rgrid, θgrid) 

        #jacobian of the local to global transformation.
        jac = dr * dθ / 4


        #computes the contribution to the W and I matrices.
        W_and_I!(W, I, met, B, prob, r, θ, ζgrid, tm)
        
        #fft the two matrices.
        p * W
        p * I

        
        #loop over the fourier components of the trial function
        for (l1, n1) in enumerate(nlist)

            #adjust the basis functions to the current coordinates/mode numbers considered.
            create_local_basis!(Φ, S, grids.θ.pf, n1, dr, dθ)

            for (l2, n2) in enumerate(nlist)

                #negatives for conjugate of test function
                create_local_basis!(Ψ, S, -grids.θ.pf, -n2, dr, dθ)

                #extract the relevant indicies from the ffted matrices.
                nind = mod(l1-l2 + Nζ, Nζ) + 1

                #loop over the Hermite elements for the trial function
                for trialr in 1:4, trialθ in 1:4
                    
                    #determines the matrix index for the trial function
                    right_ind = grid_to_index(i, j, l1, trialr, trialθ, grids)

                    #and for the test function
                    for testr in 1:4, testθ in 1:4

                        #and for the test function.
                        left_ind = grid_to_index(i, j, l2, testr, testθ, grids)

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
                                Wsum = @views gauss_integrate(Ψ[testr, testθ, :, :, :], Φ[trialr, trialθ, :, :, :], W[:, :, :, :, nind], wgr, wgθ, jac, grids.r.gp, grids.θ.gp)

                                Isum = @views gauss_integrate(Ψ[testr, testθ, :, :, :], Φ[trialr, trialθ, :, :, :], I[:, :, :, :, nind], wgr, wgθ, jac, grids.r.gp, grids.θ.gp)

                                #adds the local contribution to the global structure
                                push!(rows, left_ind)
                                push!(cols, right_ind)
                                push!(Wdata, Wsum)
                                push!(Idata, Isum)
                            end
                        else
                            
                            #integrate the local contribution to our matrices.
                            Wsum = @views gauss_integrate(Ψ[testr, testθ, :, :, :], Φ[trialr, trialθ, :, :, :], W[:, :, :, :, nind], wgr, wgθ, jac, grids.r.gp, grids.θ.gp)

                            Isum = @views gauss_integrate(Ψ[testr, testθ, :, :, :], Φ[trialr, trialθ, :, :, :], I[:, :, :, :, nind], wgr, wgθ, jac, grids.r.gp, grids.θ.gp)

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

    return Wmat, Imat
end



"""
    construct(prob::ProblemT, grids::FFFGridsT)

Constructs the two matrices using the WeakForm of the SAW governing equation. Uses Finite elements with cubic Hermite polynomials in r, θ and ζ. Returns two sparse matrices.

### Args
prob::ProblemT - Struct containing the functions and parameters that define the problem we are solving
grids::FFFGridT - Grids to solve over.
"""
function construct(prob::ProblemT, grids::FFFGridsT)

    #instantiate the grids into arrays. 
    rgrid, θgrid, ζgrid = inst_grids(grids)


    #initialise the two structs to store the metric and the magnetic field.
    met = init_empty_met()
    B = init_empty_B()

    #compute the gaussian qudrature points for finite elements.
    ξr, wgr = gausslegendre(grids.r.gp) #same as python!
    ξθ, wgθ = gausslegendre(grids.θ.gp)
    ξζ, wgζ = gausslegendre(grids.ζ.gp)

    #Gets the Hermite basis for the grids
    S = hermite_basis(ξr, ξθ, ξζ)

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

    #initialises a struct storing temporary matrices used in the weak form.
    tm = init_tm()


    #main loop, θ, ζ go to N for periodicity.
    for i in 1:grids.r.N-1, j in 1:grids.θ.N, k in 1:grids.ζ.N 


        #takes the local ξ arrays to a global arrays around the grid point.
        r, θ, ζ, dr, dθ, dζ = local_to_global(i, j, k, ξr, ξθ, ξζ, rgrid, θgrid, ζgrid) 

        #jacobian of the local to global transformation.
        jac = dr * dθ * dζ / 8 

        #computes the contribution to the W and I matrices.
        W_and_I!(W, I, met, B, prob, r, θ, ζ, tm)
        
         #adjust the basis functions to the current coordinates/mode numbers considered.
        create_local_basis!(Φ, S, grids.θ.pf, grids.ζ.pf, dr, dθ, dζ)
        #negatives for conjugate of test function
        create_local_basis!(Ψ, S, -grids.θ.pf, -grids.ζ.pf, dr, dθ, dζ)


        #loop over the Hermite elements for the trial function
        for trialr in 1:4, trialθ in 1:4, trialζ in 1:4

            #determines the matrix index for the trial function
            right_ind = grid_to_index(i, j, k, trialr, trialθ, trialζ, grids)

            #and for the test function
            for testr in 1:4, testθ in 1:4, testζ in 1:4

                 #and for the test function. 
                left_ind = grid_to_index(i, j, k, testr, testθ, testζ, grids)

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
                        Wsum = @views gauss_integrate(Ψ[testr, testθ, testζ, :, :, :, :], Φ[trialr, trialθ, trialζ, :, :, :, :], W, wgr, wgθ, wgζ, jac, grids.r.gp, grids.θ.gp, grids.ζ.gp)

                        Isum = @views gauss_integrate(Ψ[testr, testθ, testζ, :, :, :, :], Φ[trialr, trialθ, trialζ, :, :, :, :], I, wgr, wgθ, wgζ, jac, grids.r.gp, grids.θ.gp, grids.ζ.gp)


                        #adds the local contribution to the global structure
                        push!(Wdata, Wsum)
                        push!(Idata, Isum)
                        push!(rows, left_ind)
                        push!(cols, right_ind)
                    end
                else
                    
                    #integrate the local contribution to our matrices.
                    Wsum = @views gauss_integrate(Ψ[testr, testθ, testζ, :, :, :, :], Φ[trialr, trialθ, trialζ, :, :, :, :], W, wgr, wgθ, wgζ, jac, grids.r.gp, grids.θ.gp, grids.ζ.gp)

                    Isum = @views gauss_integrate(Ψ[testr, testθ, testζ, :, :, :, :], Φ[trialr, trialθ, trialζ, :, :, :, :], I, wgr, wgθ, wgζ, jac, grids.r.gp, grids.θ.gp, grids.ζ.gp)

                    #adds the local contribution to the global structure
                    push!(Wdata, Wsum)
                    push!(Idata, Isum)
                    push!(rows, left_ind)
                    push!(cols, right_ind)
                    
                end
            end
        end

    end

    #construct the sparse matrix.
    Wmat = sparse(rows, cols, Wdata)
    Imat = sparse(rows, cols, Idata)

    return Wmat, Imat
end
