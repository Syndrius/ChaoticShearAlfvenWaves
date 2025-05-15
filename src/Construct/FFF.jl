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
    met = MetT()
    B = BFieldT()

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
    tm = TM()


    #main loop, θ, ζ go to N for periodicity.
    for i in 1:grids.r.N-1, j in 1:grids.θ.N, k in 1:grids.ζ.N 

        #takes the local ξ arrays to a global arrays around the grid point.
        r, θ, ζ, dr, dθ, dζ = local_to_global(i, j, k, ξr, ξθ, ξζ, rgrid, θgrid, ζgrid) 

        #jacobian of the local to global transformation.
        jac = dr * dθ * dζ / 8 

        #computes the contribution to the W and I matrices.
        W_and_I!(W, I, B, met, prob, r, θ, ζ, tm)
        
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


"""
    construct(prob::ProblemT, grids::FFFGridsT, surfs::Array{QFMSurfaceT})

Constructs the W and I matrices using qfm surfaces to construct chaotic coordinates.
"""
function construct(prob::ProblemT, grids::FFFGridsT, surfs::Array{QFMSurfaceT})

    #vars in here really should be (s, ϑ, φ) as the grids are the new coordinates
    #however this causes issues with the data structures.

    #instantiate the grids into arrays. 
    rgrid, θgrid, ζgrid = inst_grids(grids)

    #initialise the two structs to store the metric and the magnetic field.
    tor_met = MetT()
    qfm_met = MetT()
    tor_B = BFieldT()
    qfm_B = BFieldT()

    #creates the interpolations for the surfaces.
    surf_itp, sd = create_surf_itp(surfs)

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
    tm = TM()

    #struct for storing the intermediate data for the coordinate transform
    CT = CoordTsfmT()

    #main loop, θ, ζ go to N for periodicity.
    for i in 1:grids.r.N-1, j in 1:grids.θ.N, k in 1:grids.ζ.N 


        #takes the local ξ arrays to a global arrays around the grid point.
        r, θ, ζ, dr, dθ, dζ = local_to_global(i, j, k, ξr, ξθ, ξζ, rgrid, θgrid, ζgrid) 

        #jacobian of the local to global transformation.
        jac = dr * dθ * dζ / 8 

        #computes the contribution to the W and I matrices.
        W_and_I!(W, I, tor_B, tor_met, qfm_B, qfm_met, prob, r, θ, ζ, tm, surf_itp, CT, sd)
        
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
