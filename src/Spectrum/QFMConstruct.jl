
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


function construct(prob::ProblemT, grids::FFSGridsT, surfs::Array{QFMSurfaceT})
    
    #instantiate the grids into arrays.
    rgrid, θgrid, ζgrid = inst_grids(grids)

    #for spectral method we need the length of the array
    Nζ = length(ζgrid)
    #and the mode list.
    nlist = mode_list(grids.ζ)


    #initialise the two structs to store the metric and the magnetic field.
    tor_met = MetT()
    qfm_met = MetT()
    tor_B = BFieldT()
    qfm_B = BFieldT()

    #creates the interpolations for the surfaces.
    surf_itp = create_surf_itp(surfs)

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

    #struct for storing the intermediate data for the coordinate transform
    CT = CoordTsfmT()

    
    #main loop, θ goes to N for periodicity.
    for i in 1:grids.r.N-1, j in 1:grids.θ.N 

        #takes the local ξ arrays to a global arrays around the grid point.
        r, θ, dr, dθ = local_to_global(i, j, ξr, ξθ, rgrid, θgrid) 

        #jacobian of the local to global transformation.
        jac = dr * dθ / 4


        #computes the contribution to the W and I matrices.
        #W_and_I!(W, I, met, B, prob, r, θ, ζgrid, tm)
        W_and_I!(W, I, tor_met, tor_B, qfm_met, qfm_B, prob, r, θ, ζgrid, tm, surf_itp, CT)
        
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




#attempt at continuum version, structure of evals looks okish, but magnitude is wild, and seems to scale with M, N.
function qfm_continuum(prob::ProblemT, grids::ContGridsT, surfs::Array{QFMSurfaceT})

    #does the continuum need to use fourier spectral??

    #probably not, but that is what we do.

    #so this is actually returning normal values.
    #unsure how the size of this will correlate with the surface size
    #and/or output size.
    #i.e. does the θ and ζ size have to match the surface size?
    #probably not, I think the surface just gives us an s, ϑ, φ value
    #for each point in our grid?
    #may be odd how it combines with fourier spectral method tho.
    rgrid, θgrid, ζgrid = inst_grids(grids)

    Nθ = length(θgrid)
    Nζ = length(ζgrid)

    #custom grid that will work for the two wave case.
    #rgrid = LinRange(1.1095, 1.112, 10)

    #see this seems a bit odd tbh.
    #does seem like finite elements would be better here?
    #that takes fkn ages though so we would need to implement a parallel construct for continuum case... big rip.
    mlist = mode_list(grids.θ)
    nlist = mode_list(grids.ζ)

    #allows us to create the surfaces eternally, which will be useful
    #to compare how many we need.
    #this is probably a good idea long term
    #if isnone(surf.rcos_itp)
        #then create surfaces

        #thinking this struct will store all the random af params Zhisong uses.
    #   action(QFMSurfaceT)
    #end

    #for now, we will assume the surface obj is created, but we need to create the interpolation still.
    surf_itp = create_surf_itp(surfs);

    tor_met = MetT()
    qfm_met = MetT()
    tor_B = BFieldT()
    qfm_B = BFieldT()

    mat_size = matrix_size(grids)

    ωlist = zeros(grids.r.N, mat_size)

    #pretty sure we have functions to do a bunch of this shite.
    I = zeros(ComplexF64, 9, 9, 1, Nθ, Nζ)
    W = zeros(ComplexF64, 9, 9, 1, Nθ, Nζ)

    Imat = zeros(ComplexF64, mat_size, mat_size)
    Wmat = zeros(ComplexF64, mat_size, mat_size)

    #so turns out our original continuum funciton is garbage.
    Icont = zeros(1, 1, 1, Nθ, Nζ)

    Wcont = zeros(3, 3, 1, Nθ, Nζ)

    #probably not how this works.
    pI = plan_fft(Icont, [4, 5])
    pW = plan_fft(Wcont, [4, 5])

    #structs storing temporary arrays.
    CT = CoordTsfmT()
    tm = TM()

    for (i, r) in enumerate(rgrid)

        tor_met = MetT()
        qfm_met = MetT()
        tor_B = BFieldT()
        qfm_B = BFieldT()
        CT = CoordTsfmT()
        tm = TM()
        Imat = zeros(ComplexF64, mat_size, mat_size)
        Wmat = zeros(ComplexF64, mat_size, mat_size)

        W_and_I!(W, I, tor_met, tor_B, qfm_met, qfm_B, prob, [r], θgrid, ζgrid, tm, surf_itp, CT)

        #display(CT.JM)
        #display(tor_B.B)
        #display(W)

        Icont = pI * I[[1], [1], :, :, :]
        Wcont = pW * W[[1, 5, 6], [1, 5, 6], :, :, :]
        

        Imat .= 0.0
        Wmat .= 0.0

        for (k1,m1) in enumerate(mlist), (l1, n1) in enumerate(nlist)
            #left_ind = cont_grid_to_index(k1, l1, grids.ζ.count)
            left_ind = grid_to_index(k1, l1, grids)

            #equivalent to create local basis. Much simpler in contninuum case.
            Φ = [1, m1*1im, n1*1im]
            
            for (k2, m2) in enumerate(mlist), (l2, n2) in enumerate(nlist)
                right_ind = grid_to_index(k2, l2, grids)
                
                Ψ = [1, -m2*1im, -n2*1im]

                #mind = mod(m1-m2 + nθ, nθ) + 1
                #nind = mod(n1-n2 + nζ, nζ) + 1
                mind = mod(k1-k2 + Nθ, Nθ) + 1
                nind = mod(l1-l2 + Nζ, Nζ) + 1

                
                

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


        #odds that this will be Hermitian?
        #so all the matrix elements are freekin enourmous.
        #explains the large evals.
        #display(Wmat)
        #Wmat especially seems to be wopping.
        #slab vs cylinder is not the problemo!
        
        #vals = eigvals(Hermitian(Wmat), Hermitian(Imat))
        #basic benchmark test shows that this is Hermitian!
        vals = eigvals(Wmat, Imat)


        ωlist[i, :] = prob.geo.R0 * sqrt.(abs.(vals))
    end

    #should probably return an evals Obj.
    return ωlist


end
