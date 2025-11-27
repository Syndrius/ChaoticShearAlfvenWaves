"""
    construct(prob::ProblemT, grids::FFFGridsT)

Constructs the two matrices using the PeakForm of the SAP governing equation. Uses Finite elements with cubic Hermite polynomials in r, x2 and x3. Returns two sparse matrices.

### Args
prob::ProblemT - Struct containing the functions and parameters that define the problem we are solving
grids::FFFGridT - Grids to solve over.
"""
function construct(prob::ProblemT, grids::FFFGridsT)

    #instantiate the grids into arrays. 
    x1grid, x2grid, x3grid = inst_grids(grids)

    #initialise the two structs to store the metric and the magnetic field.
    met = MetT()
    B = BFieldT()

    #compute the gaussian qudrature points for finite elements.
    ξx1, ξx2, ξx3, wgx1, wgx2, wgx3 = gauss_points(grids)

    #Gets the Hermite basis for the grids
    S = hermite_basis(ξx1, ξx2, ξx3)

    #array for storing the scaling of the tangent basis functions when transforming
    #from the local ξ∈[-1, 1] to global x∈[x_i, x_{i+1}] domain
    ts = ones(size(S.H))

    #creates the trial and test function arrays.
    #these store the basis functions for each derivative
    #and finite elements basis 
    Φ = init_trial_function(grids)
    Ψ = init_trial_function(grids)

    #arrays to store the row, column and data of each matrix element
    #used for constructing sparse matrices.
    rows = Array{Int64}(undef, 0) 
    cols = Array{Int64}(undef, 0) 
    Qdata = Array{ComplexF64}(undef, 0)
    Pdata = Array{ComplexF64}(undef, 0)

    #determines the indicies in the matrices corresponding to the dirichlet boundary conditions for x1.
    boundary_inds = compute_boundary_inds(grids)

    #generalised eval problem PΦ = ω^2 Q Φ
    #these matrices store the local contribution, i.e. at each grid point, for the global matrices Q and P.
    Q = init_local_matrix(grids)
    P = init_local_matrix(grids)

    #initialises a struct storing temporary matrices used in the weak form.
    tm = TM()

    #will want a function for this!
    x1 = zeros(length(ξx1))
    x2 = zeros(length(ξx2))
    x3 = zeros(length(ξx3))
    dx = zeros(3)

    #main loop, x2, x3 go to N for periodicity.
    for i in 1:grids.x1.N-1, j in 1:grids.x2.N, k in 1:grids.x3.N 

        #takes the local ξ arrays to a global arrays around the grid point.
        #number of args probably tells us that this could be improved
        jac = local_to_global!(x1, x2, x3, dx, i, j, k, ξx1, ξx2, ξx3, x1grid, x2grid, x3grid)

        #computes the contribution to the P and Q matrices.
        weak_form!(P, Q, B, met, prob, x1, x2, x3, tm)

        #transforms the local basis function to the global.
        update_trial_function!(Φ, S, grids.x2.pf, grids.x3.pf, dx, ts)
        #negatives for conjugate of test function
        update_trial_function!(Ψ, S, -grids.x2.pf, -grids.x3.pf, dx, ts)

        #loop over the Hermite elements for the trial function
        for trialx1 in 1:4, trialx2 in 1:4, trialx3 in 1:4

            #determines the matrix index for the trial function
            right_ind = grid_to_index(i, j, k, trialx1, trialx2, trialx3, grids)

            #and for the test function
            for testx1 in 1:4, testx2 in 1:4, testx3 in 1:4

                #and for the test function. 
                left_ind = grid_to_index(i, j, k, testx1, testx2, testx3, grids)

                #only check for boundaries if this is true
                #no other i's can possibly give boundaries
                if i==1 || i==grids.x1.N-1

                    if left_ind == right_ind && left_ind in boundary_inds

                        #diagonals for boundary conditions are set to 1.
                        push!(rows, left_ind)
                        push!(cols, right_ind)
                        push!(Pdata, 1.0 + 0.0im)
                        push!(Qdata, 1.0 + 0.0im)
                        
                    #otherwise the boundaries are set to zero, which for sparse matrices
                    #is the same as leaving blank.
                    elseif left_ind in boundary_inds
                        continue
                    elseif right_ind in boundary_inds
                        continue
                    #otherwise a regular case for these indicies.
                    else
                        
                        #integrate the local contribution to our matrices.
                        Psum = @views gauss_integrate(Ψ[testx1, testx2, testx3, :, :, :, :], Φ[trialx1, trialx2, trialx3, :, :, :, :], P, wgx1, wgx2, wgx3, jac, grids.x1.gp, grids.x2.gp, grids.x3.gp)

                        Qsum = @views gauss_integrate(Ψ[testx1, testx2, testx3, :, :, :, :], Φ[trialx1, trialx2, trialx3, :, :, :, :], Q, wgx1, wgx2, wgx3, jac, grids.x1.gp, grids.x2.gp, grids.x3.gp)

                        #adds the local contribution to the global structure
                        push!(Pdata, Psum)
                        push!(Qdata, Qsum)
                        push!(rows, left_ind)
                        push!(cols, right_ind)
                    end
                else
                    
                    #integrate the local contribution to our matrices.
                    Psum = @views gauss_integrate(Ψ[testx1, testx2, testx3, :, :, :, :], Φ[trialx1, trialx2, trialx3, :, :, :, :], P, wgx1, wgx2, wgx3, jac, grids.x1.gp, grids.x2.gp, grids.x3.gp)

                    Qsum = @views gauss_integrate(Ψ[testx1, testx2, testx3, :, :, :, :], Φ[trialx1, trialx2, trialx3, :, :, :, :], Q, wgx1, wgx2, wgx3, jac, grids.x1.gp, grids.x2.gp, grids.x3.gp)

                    #adds the local contribution to the global structure
                    push!(Pdata, Psum)
                    push!(Qdata, Qsum)
                    push!(rows, left_ind)
                    push!(cols, right_ind)

                end
            end
        end

    end

    #construct the sparse matrix.
    Pmat = sparse(rows, cols, Pdata)
    Qmat = sparse(rows, cols, Qdata)

    return Pmat, Qmat
end

"""
    construct(prob::ProblemT, grids::FFFGridsT, surfs::Array{QFMSurfaceT})

Constructs the P and Q matrices using qfm surfaces to construct chaotic coordinates.
"""
function construct(prob::ProblemT, grids::FFFGridsT, surfs::Array{QFMSurfaceT})

    #instantiate the grids into arrays. 
    x1grid, x2grid, x3grid = inst_grids(grids)

    #initialise the two structs to store the metric and the magnetic field.
    #one set for each coordinate system
    tor_met = MetT()
    qfm_met = MetT()
    tor_B = BFieldT()
    qfm_B = BFieldT()

    #creates the interpolations for the surfaces.
    surf_itp, sd = create_surf_itp(surfs)

    #compute the gaussian qudrature points for finite elements.
    ξx1, ξx2, ξx3, wgx1, wgx2, wgx3 = gauss_points(grids)

    #Gets the Hermite basis for the grids
    S = hermite_basis(ξx1, ξx2, ξx3)

    #array for storing the scaling of the tangent basis functions when transforming
    #from the local ξ∈[-1, 1] to global x∈[x_i, x_{i+1}] domain
    ts = ones(size(S.H))

    #creates the trial and test function arrays.
    #these store the basis functions for each derivative
    #and finite elements basis 
    Φ = init_trial_function(grids)
    Ψ = init_trial_function(grids)

    #arrays to store the row, column and data of each matrix element
    #used for constructing sparse matrices.
    rows = Array{Int64}(undef, 0) 
    cols = Array{Int64}(undef, 0) 
    Qdata = Array{ComplexF64}(undef, 0)
    Pdata = Array{ComplexF64}(undef, 0)

    #determines the indicies in the matrices corresponding to the dirichlet boundary conditions for x1.
    boundary_inds = compute_boundary_inds(grids)

    #generalised eval problem PΦ = ω^2 Q Φ
    #these matrices store the local contribution, i.e. at each grid point, for the global matrices Q and P.
    Q = init_local_matrix(grids)
    P = init_local_matrix(grids)

    #initialises a struct storing temporary matrices used in the weak form.
    tm = TM()

    #will want a function for this!
    x1 = zeros(length(ξx1))
    x2 = zeros(length(ξx2))
    x3 = zeros(length(ξx3))
    dx = zeros(3)

    #struct for storing the intermediate data for the coordinate transform
    CT = CoordTransformT()

    #main loop, x2, x3 go to N for periodicity.
    for i in 1:grids.x1.N-1, j in 1:grids.x2.N, k in 1:grids.x3.N 

        #takes the local ξ arrays to a global arrays around the grid point.
        #number of args probably tells us that this could be improved
        jac = local_to_global!(x1, x2, x3, dx, i, j, k, ξx1, ξx2, ξx3, x1grid, x2grid, x3grid)

        #computes the contribution to the P and Q matrices.
        weak_form!(P, Q, tor_B, tor_met, qfm_B, qfm_met, prob, x1, x2, x3, tm, surf_itp, CT, sd)

        #transforms the local basis function to the global.
        update_trial_function!(Φ, S, grids.x2.pf, grids.x3.pf, dx, ts)
        #negatives for conjugate of test function
        update_trial_function!(Ψ, S, -grids.x2.pf, -grids.x3.pf, dx, ts)

        #loop over the Hermite elements for the trial function
        for trialx1 in 1:4, trialx2 in 1:4, trialx3 in 1:4

            #determines the matrix index for the trial function
            right_ind = grid_to_index(i, j, k, trialx1, trialx2, trialx3, grids)

            #and for the test function
            for testx1 in 1:4, testx2 in 1:4, testx3 in 1:4

                #and for the test function. 
                left_ind = grid_to_index(i, j, k, testx1, testx2, testx3, grids)

                #only check for boundaries if this is true
                #no other i's can possibly give boundaries
                if i==1 || i==grids.x1.N-1

                    if left_ind == right_ind && left_ind in boundary_inds

                        #diagonals for boundary conditions are set to 1.
                        push!(rows, left_ind)
                        push!(cols, right_ind)
                        push!(Pdata, 1.0 + 0.0im)
                        push!(Qdata, 1.0 + 0.0im)
                        
                    #otherwise the boundaries are set to zero, which for sparse matrices
                    #is the same as leaving blank.
                    elseif left_ind in boundary_inds
                        continue
                    elseif right_ind in boundary_inds
                        continue
                    #otherwise a regular case for these indicies.
                    else
                        
                        #integrate the local contribution to our matrices.
                        Psum = @views gauss_integrate(Ψ[testx1, testx2, testx3, :, :, :, :], Φ[trialx1, trialx2, trialx3, :, :, :, :], P, wgx1, wgx2, wgx3, jac, grids.x1.gp, grids.x2.gp, grids.x3.gp)

                        Qsum = @views gauss_integrate(Ψ[testx1, testx2, testx3, :, :, :, :], Φ[trialx1, trialx2, trialx3, :, :, :, :], Q, wgx1, wgx2, wgx3, jac, grids.x1.gp, grids.x2.gp, grids.x3.gp)

                        #adds the local contribution to the global structure
                        push!(Pdata, Psum)
                        push!(Qdata, Qsum)
                        push!(rows, left_ind)
                        push!(cols, right_ind)
                    end
                else
                    
                    #integrate the local contribution to our matrices.
                    Psum = @views gauss_integrate(Ψ[testx1, testx2, testx3, :, :, :, :], Φ[trialx1, trialx2, trialx3, :, :, :, :], P, wgx1, wgx2, wgx3, jac, grids.x1.gp, grids.x2.gp, grids.x3.gp)

                    Qsum = @views gauss_integrate(Ψ[testx1, testx2, testx3, :, :, :, :], Φ[trialx1, trialx2, trialx3, :, :, :, :], Q, wgx1, wgx2, wgx3, jac, grids.x1.gp, grids.x2.gp, grids.x3.gp)

                    #adds the local contribution to the global structure
                    push!(Pdata, Psum)
                    push!(Qdata, Qsum)
                    push!(rows, left_ind)
                    push!(cols, right_ind)

                end
            end
        end

    end

    #construct the sparse matrix.
    Pmat = sparse(rows, cols, Pdata)
    Qmat = sparse(rows, cols, Qdata)

    return Pmat, Qmat
end

