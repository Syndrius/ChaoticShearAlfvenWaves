"""
    construct(prob::ProblemT, grids::FSSGridsT)

Constructs the two matrices using the WeakForm of the SAW governing equation. Uses Finite elements with cubic Hermite polynomials in r, and the fourier spectral method in θ and ζ. Returns two sparse matrices.

### Args
prob::ProblemT - Struct containing the functions and parameters that define the problem we are solving
grids::FSSGridT - Grids to solve over.
"""
function construct(prob::ProblemT, grids::FSSGridsT)

    #removing the kwargs from here seems to have helped with the weird warnings, I think they are just a symptom of a previously compiled version!
    #island being nothing is stupid, either need to separate cases with and without island, or construct an empty (A=0) island if it is nothing.

    #nθ, mlist, θgrid = spectral_grid(grids.pmd)
    #nζ, nlist, ζgrid = spectral_grid(grids.tmd)

    rgrid, θgrid, ζgrid = inst_grids(grids)

    Nθ = length(θgrid)
    Nζ = length(ζgrid)

    mlist = mode_list(grids.θ)
    nlist = mode_list(grids.ζ)



    #initialise the two structs.
    met = init_empty_met()
    B = init_empty_B()

    ξ, wg = gausslegendre(grids.r.gp) #same as python!

    #gets the basis 
    H, dH, ddH = hermite_basis(ξ)

    #the trial function
    #4 is the number of Hermite shape functions
    #10 is Φ and all its relevant derivatives.
    #under new method we don't use the zeroth derivative, so these could be replaced with 9 
    #TODO shape of this needs a proper description!

    #probably want a function for this!
    Φ = init_bases_function(grids)
    #the test function.
    Ψ = init_bases_function(grids)
    


    #generalised eval problem WΦ = ω^2 I Φ
    #Imat = zeros(ComplexF64, 2 * rd.grid_size * pmd.count * tmd.count, 2 * rd.grid_size * pmd.count * tmd.count)
    #Wmat = zeros(ComplexF64, 2 * rd.grid_size * pmd.count * tmd.count, 2 * rd.grid_size * pmd.count * tmd.count)

    #probably possible to know the size of this first.
    #can probably determine the maximum sized int we need based on size of matrix.
    #need to determine the size of these bad bois, they are eating up a bit of time.
    #this seems to have made a minimal difference, but is surely better practise right.
    #arr_length = compute_length(grids.r.N, grids.θ.count, grids.ζ.count)

    #arr_count = 1 #this may be wrong... gives error for v small matrix...

    #rows = Array{Int64}(undef, arr_length) #ints
    #cols = Array{Int64}(undef, arr_length) #ints
    #Idata = Array{ComplexF64}(undef, arr_length)
    #Wdata = Array{ComplexF64}(undef, arr_length)

    rows = Array{Int64}(undef, 0) #ints
    cols = Array{Int64}(undef, 0) #ints
    Idata = Array{ComplexF64}(undef, 0)
    Wdata = Array{ComplexF64}(undef, 0)


    #either need a condition in case m=0 or just an error message.
    #boundary_inds = compute_boundary_inds(grids.r.N, grids.θ.count, grids.ζ.count, collect(mlist))
    boundary_inds = compute_boundary_inds(grids)
    #display(size(boundary_inds))
    #display(boundary_inds)

    #feking awful name
    #not sure how to desribe what I and W here actually are.
    I = local_matrix_size(grids)
    W = local_matrix_size(grids)


    #this gives a warning but seems to work perfectly
    #and incredibly efficiently!
    #this essentially defines a plane to fft I, (and W as they are the same size), which can be exectued
    #via p * I, this is done in place and seems to be mad efficient.
    p = plan_fft!(W, [4, 5])

    #Isum = 0.0 + 0.0im
    #Wsum = 0.0 + 0.0im

    #bounds_count = 0

    #now we loop through the grid

    #rgrid = construct_rgrid(grids)

    for i in 1:grids.r.N-1

        r, dr = local_to_global(i, ξ, rgrid)

        jac = dr/2 #following thesis!


        #I_and_W!(I, W, B, q_profile, met, compute_met, dens, r, θgrid, ζgrid, δ, isl, R0)

        W_and_I!(W, I, met, B, prob, r, θgrid, ζgrid)
        #W_tor, I_tor = stupid_W_and_I!(W, I, met, B, prob, r, θgrid, ζgrid)
        #stupid_W_and_I!(W, I, met, B, prob, r, θgrid, ζgrid)

        #display(W[:, :, 1, 1, 1])


        #uses the fft plan to take the fft of our two matrices.
        p * W
        p * I
        #p * W_tor
        #p * I_tor

        #W[1, 1, :, :, :] .= W_tor[1, 1, :, :, :]
        #I[1, 1, :, :, :] .= I_tor[1, 1, :, :, :]
        #I[4, :, :, :, :] .= I_tor[4, :, :, :, :]
        #W[4, :, :, :, :] .= W_tor[4, :, :, :, :]
        #I[:, 4, :, :, :] .= I_tor[:, 4, :, :, :]
        #W[:, 4, :, :, :] .= W_tor[:, 4, :, :, :]

        
        #loop over the fourier components of the trial function
        for (k1,m1) in enumerate(mlist), (l1, n1) in enumerate(nlist)

            create_local_basis!(Φ, H, dH, ddH, m1, n1, jac)

            for (k2, m2) in enumerate(mlist), (l2, n2) in enumerate(nlist)

                #negatives for conjugate
                create_local_basis!(Ψ, H, dH, ddH, -m2, -n2, jac)

                #extract the relevant indicies from the ffted matrices.
                mind = mod(k1-k2 + Nθ, Nθ) + 1
                nind = mod(l1-l2 + Nζ, Nζ) + 1


                for trialsf in 1:4

                    right_ind = grid_to_index(i, k1, l1, trialsf, grids)

                    for testsf in 1:4
                        #display("testsf")
                        #display(testsf)

                        
                        left_ind = grid_to_index(i, k2, l2, testsf, grids)

                        #only check for boundaries if this is true
                        #no other i's can possibly give boundaries
                        if i==1 || i==grids.r.N-1


                            if left_ind == right_ind && left_ind in boundary_inds

                                #rows[arr_count] = left_ind
                                #cols[arr_count] = right_ind
                                #Wdata[arr_count] = 1.0 + 0.0im
                                #Idata[arr_count] = 1.0 + 0.0im

                                push!(rows, left_ind)
                                push!(cols, right_ind)

                                push!(Wdata, 1.0 + 0.0im)
                                push!(Idata, 1.0 + 0.0im)
                                
                                #arr_count += 1

                                #bounds_count += 1
                            
                            #otherwise the boundaries are set to zero, which for sparse matrices
                            #is the same as leaving blank.
                            elseif left_ind in boundary_inds
                                continue
                            elseif right_ind in boundary_inds
                                continue
                            #otherwise a regular case for these indicies.
                            else
                                #rows[arr_count] = left_ind
                                #cols[arr_count] = right_ind
                                push!(rows, left_ind)
                                push!(cols, right_ind)

                                
                                

                                Wsum = @views gauss_integrate(Ψ[testsf, :, :], Φ[trialsf, :, :], W[:, :, :, mind, nind], wg, jac, grids.r.gp)


                                Isum = @views gauss_integrate(Ψ[testsf, :, :], Φ[trialsf, :, :], I[:, :, :, mind, nind], wg, jac, grids.r.gp)

                                #Wdata[arr_count] = Wsum
                                #Idata[arr_count] = Isum
                                push!(Wdata, Wsum)
                                push!(Idata, Isum)
                                
                                #arr_count += 1
                            end
                        else
                            
                            #rows[arr_count] = left_ind
                            #cols[arr_count] = right_ind

                            push!(rows, left_ind)
                            push!(cols, right_ind) 
                                

                            Wsum = @views gauss_integrate(Ψ[testsf, :, :], Φ[trialsf, :, :], W[:, :, :, mind, nind], wg, jac, grids.r.gp)


                            Isum = @views gauss_integrate(Ψ[testsf, :, :], Φ[trialsf, :, :], I[:, :, :, mind, nind], wg, jac, grids.r.gp)

                            #Wdata[arr_count] = Wsum
                            #Idata[arr_count] = Isum

                            push!(Wdata, Wsum)
                            push!(Idata, Isum)
                            
                            #arr_count += 1
                        end
                    end
                end

            end

        end
    end

    #maybe more consisnt for this function to return the rows and data as per parallal case.
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


    #island being nothing is stupid, either need to separate cases with and without island, or construct an empty (A=0) island if it is nothing.

    

    rgrid, θgrid, ζgrid = inst_grids(grids)

    Nζ = length(ζgrid)
    nlist = mode_list(grids.ζ)


    #initialise the two structs.
    met = init_empty_met()
    #MetT(zeros(3, 3), zeros(3, 3), zeros(3, 3, 3), zeros(3, 3, 3), 0.0, zeros(3))
    B = init_empty_B()
    #FieldT(zeros(3), zeros(3), zeros(3, 3), zeros(3, 3), 0.0, zeros(3))

    #not sure if this should be combined into 1 or something, focus on getting to work first.
    ξr, wgr = gausslegendre(grids.r.gp) #same as python!
    ξθ, wgθ = gausslegendre(grids.θ.gp)

    #gets the basis 

    S, dSr, dSθ, ddSrr, ddSrθ, ddSθθ = hermite_basis(ξr, ξθ)

    #the trial function
    #4 is the number of Hermite shape functions
    #10 is Φ and all its relevant derivatives.
    #under new method we don't use the zeroth derivative, so these could be replaced with 9 

    #shape of this will be cooked, expect 4-> 16, unsure if we combined rd and θd yet, leave separate for now.
    Φ = init_bases_function(grids)
    #the test function.
    Ψ = init_bases_function(grids)
    


    #generalised eval problem WΦ = ω^2 I Φ
    #Imat = zeros(ComplexF64, 2 * rd.grid_size * pmd.count * tmd.count, 2 * rd.grid_size * pmd.count * tmd.count)
    #Wmat = zeros(ComplexF64, 2 * rd.grid_size * pmd.count * tmd.count, 2 * rd.grid_size * pmd.count * tmd.count)

    #probably possible to know the size of this first.
    #can probably determine the maximum sized int we need based on size of matrix.
    #need to determine the size of these bad bois, they are eating up a bit of time.
    #this seems to have made a minimal difference, but is surely better practise right.
    #arr_length = compute_length(grids.rd.N, grids.pmd.count, grids.tmd.count)

    #arr_count = 1 #this may be wrong... gives error for v small matrix...

    #probably won't know the lengths anymore!
    #rows = Array{Int64}(undef, arr_length) #ints
    #cols = Array{Int64}(undef, arr_length) #ints
    #Idata = Array{ComplexF64}(undef, arr_length)
    #Wdata = Array{ComplexF64}(undef, arr_length)

    rows = Array{Int64}(undef, 0) #ints
    cols = Array{Int64}(undef, 0) #ints
    Idata = Array{ComplexF64}(undef, 0)
    Wdata = Array{ComplexF64}(undef, 0)


    #either need a condition in case m=0 or just an error message.
    boundary_inds = compute_boundary_inds(grids)
    #display(size(boundary_inds))
    #display(boundary_inds)

    #these will hopefully be smaller I think!
    I = local_matrix_size(grids)
    W = local_matrix_size(grids)


    #this gives a warning but seems to work perfectly
    #and incredibly efficiently!
    #this essentially defines a plane to fft I, (and W as they are the same size), which can be exectued
    #via p * I, this is done in place and seems to be mad efficient.
    #only a single ft now!
    p = plan_fft!(W, [5])

    #Isum = 0.0 + 0.0im
    #Wsum = 0.0 + 0.0im


    #now we loop through the grid

    #will we want a clustered θgrid??? probably not???
    #but we can generalise this function as well if we want it to work with θ, even without clustering!
    #rgrid, θgrid = construct_grids_zf(grids)

    for i in 1:grids.r.N-1, j in 1:grids.θ.N #go to N for periodicity!


        r, θ, dr, dθ = local_to_global(i, j, ξr, ξθ, rgrid, θgrid) #wot is θgrid? will need to be constructed I think.

        #Hopefully this is correct!
        jac = dr * dθ / 4 #following thesis!


        #I_and_W!(I, W, B, q_profile, met, compute_met, dens, r, θgrid, ζgrid, δ, isl, R0)

        #hopefully this step will be smaller! but we have twice the loop, so everything else will be longer!
        W_and_I!(W, I, met, B, prob, r, θ, ζgrid)
        #W_tor, I_tor = stupid_W_and_I!(W, I, met, B, prob, r, θgrid, ζgrid)
        #stupid_W_and_I!(W, I, met, B, prob, r, θgrid, ζgrid)

        #display(W[:, :, 1, 1, 1])


        #uses the fft plan to take the fft of our two matrices.
        p * W
        p * I


        
        #loop over the fourier components of the trial function
        for (l1, n1) in enumerate(nlist)

            #note we haven't implemented pf for r, seems pointless.
            create_local_basis!(Φ, S, dSr, dSθ, ddSrr, ddSrθ, ddSθθ, grids.θ.pf, n1, dr, dθ)

            for (l2, n2) in enumerate(nlist)

                #negatives for conjugate, will assume the phase factor is conjugate as well.
                create_local_basis!(Ψ, S, dSr, dSθ, ddSrr, ddSrθ, ddSθθ, -grids.θ.pf, -n2, dr, dθ)

                #extract the relevant indicies from the ffted matrices.
                #mind = mod(k1-k2 + nθ, nθ) + 1
                nind = mod(l1-l2 + Nζ, Nζ) + 1


                for trialr in 1:4, trialθ in 1:4
                    
                    #may need a θN or something!
                    right_ind = grid_to_index(i, j, l1, trialr, trialθ, grids)

                    for testr in 1:4, testθ in 1:4
                        #display("testsf")
                        #display(testsf)

                        
                        left_ind = grid_to_index(i, j, l2, testr, testθ, grids)

                        #only check for boundaries if this is true
                        #no other i's can possibly give boundaries
                        
                        if i==1 || i==grids.r.N-1


                            if left_ind == right_ind && left_ind in boundary_inds

                                #rows[arr_count] = left_ind
                                #cols[arr_count] = right_ind
                                #Wdata[arr_count] = 1.0 + 0.0im
                                #Idata[arr_count] = 1.0 + 0.0im

                                push!(rows, left_ind)
                                push!(cols, right_ind)
                                #TODO not sure about this!
                                push!(Wdata, 0.5 + 0.0im)
                                push!(Idata, 0.5 + 0.0im)

                                #display((i, j))
                                
                                #arr_count += 1

                                #bounds_count += 1
                            
                            #otherwise the boundaries are set to zero, which for sparse matrices
                            #is the same as leaving blank.
                            elseif left_ind in boundary_inds
                                continue
                            elseif right_ind in boundary_inds
                                continue
                            #otherwise a regular case for these indicies.
                            else
                                #rows[arr_count] = left_ind
                                #cols[arr_count] = right_ind
                                push!(rows, left_ind)
                                push!(cols, right_ind)
                                


                                Wsum = @views gauss_integrate(Ψ[testr, testθ, :, :, :], Φ[trialr, trialθ, :, :, :], W[:, :, :, :, nind], wgr, wgθ, jac, grids.r.gp, grids.θ.gp)


                                Isum = @views gauss_integrate(Ψ[testr, testθ, :, :, :], Φ[trialr, trialθ, :, :, :], I[:, :, :, :, nind], wgr, wgθ, jac, grids.r.gp, grids.θ.gp)

                                push!(Wdata, Wsum)
                                push!(Idata, Isum)
                            end
                        else
                            
                            #rows[arr_count] = left_ind
                            #cols[arr_count] = right_ind
                            push!(rows, left_ind)
                            push!(cols, right_ind)
                                

                            Wsum = @views gauss_integrate(Ψ[testr, testθ, :, :, :], Φ[trialr, trialθ, :, :, :], W[:, :, :, :, nind], wgr, wgθ, jac, grids.r.gp, grids.θ.gp)


                            Isum = @views gauss_integrate(Ψ[testr, testθ, :, :, :], Φ[trialr, trialθ, :, :, :], I[:, :, :, :, nind], wgr, wgθ, jac, grids.r.gp, grids.θ.gp)

                            push!(Wdata, Wsum)
                            push!(Idata, Isum)
                            
                        end
                    end
                end

            end

        end
    end


    #maybe more consisnt for this function to return the rows and data as per parallal case.
    Wmat = sparse(rows, cols, Wdata)
    Imat = sparse(rows, cols, Idata)

    #display(Matrix(Wmat))

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

    #island being nothing is stupid, either need to separate cases with and without island, or construct an empty (A=0) island if it is nothing.

    
    rgrid, θgrid, ζgrid = inst_grids(grids)


    #initialise the two structs.
    met = init_empty_met()
    B = init_empty_B()

    #not sure if this should be combined into 1 or something, focus on getting to work first.
    ξr, wgr = gausslegendre(grids.r.gp) #same as python!
    ξθ, wgθ = gausslegendre(grids.θ.gp)
    ξζ, wgζ = gausslegendre(grids.ζ.gp)

    #gets the basis 
    #S = hermite_basis(ξr, ξθ, ξζ)
    #ideally these would be combined in some way, this is fkn stupid.
    S, dSr, dSθ, dSζ, ddSrr, ddSrθ, ddSrζ, ddSθθ, ddSθζ, ddSζζ = hermite_basis(ξr, ξθ, ξζ)

    #the trial function
    #4 is the number of Hermite shape functions
    #10 is Φ and all its relevant derivatives.
    #under new method we don't use the zeroth derivative, so these could be replaced with 9 

    #order of this is extemely important for @views.
    #because we integrate over each basis individually, the basis indicies (the 4's), should go first,
    #then when we use @views, the entire block of memory is combined for more efficient summation.
    #probably need a proper comment description for this structure as it is cooked beyond belief.
    Φ = init_bases_function(grids)
    #the test function.
    Ψ = init_bases_function(grids)


    #generalised eval problem WΦ = ω^2 I Φ
    #Imat = zeros(ComplexF64, 2 * rd.grid_size * pmd.count * tmd.count, 2 * rd.grid_size * pmd.count * tmd.count)
    #Wmat = zeros(ComplexF64, 2 * rd.grid_size * pmd.count * tmd.count, 2 * rd.grid_size * pmd.count * tmd.count)

    #probably possible to know the size of this first.
    #can probably determine the maximum sized int we need based on size of matrix.
    #need to determine the size of these bad bois, they are eating up a bit of time.
    #this seems to have made a minimal difference, but is surely better practise right.
    #arr_length = compute_length(grids.rd.N, grids.pmd.count, grids.tmd.count)

    #arr_count = 1 #this may be wrong... gives error for v small matrix...

    #probably won't know the lengths anymore!
    #rows = Array{Int64}(undef, arr_length) #ints
    #cols = Array{Int64}(undef, arr_length) #ints
    #Idata = Array{ComplexF64}(undef, arr_length)
    #Wdata = Array{ComplexF64}(undef, arr_length)

    rows = Array{Int64}(undef, 0) #ints
    cols = Array{Int64}(undef, 0) #ints
    Idata = Array{ComplexF64}(undef, 0)
    Wdata = Array{ComplexF64}(undef, 0)


    #either need a condition in case m=0 or just an error message.
    boundary_inds = compute_boundary_inds(grids)
    #display(size(boundary_inds))
    #display(boundary_inds)

    #these will hopefully be smaller I think!
    I = local_matrix_size(grids)
    W = local_matrix_size(grids)


    #this gives a warning but seems to work perfectly
    #and incredibly efficiently!
    #this essentially defines a plane to fft I, (and W as they are the same size), which can be exectued
    #via p * I, this is done in place and seems to be mad efficient.
    #only a single ft now!


    #Isum = 0.0 + 0.0im
    #Wsum = 0.0 + 0.0im


    #now we loop through the grid

    #will we want a clustered θgrid??? probably not???
    #but we can generalise this function as well if we want it to work with θ, even without clustering!
    #rgrid, θgrid = construct_grids_zf(grids)

    for i in 1:grids.r.N-1, j in 1:grids.θ.N, k in 1:grids.ζ.N #go to N for periodicity!


        #display((i, j, k))
        r, θ, ζ, dr, dθ, dζ = local_to_global(i, j, k, ξr, ξθ, ξζ, rgrid, θgrid, ζgrid) #wot is θgrid? will need to be constructed I think.

        #Hopefully this is correct!
        jac = dr * dθ * dζ / 8 #following thesis!


        #I_and_W!(I, W, B, q_profile, met, compute_met, dens, r, θgrid, ζgrid, δ, isl, R0)

        #hopefully this step will be smaller! but we have twice the loop, so everything else will be longer!
        W_and_I!(W, I, met, B, prob, r, θ, ζ)
        #W_tor, I_tor = stupid_W_and_I!(W, I, met, B, prob, r, θgrid, ζgrid)
        #stupid_W_and_I!(W, I, met, B, prob, r, θgrid, ζgrid)

        #display(W[:, :, 1, 1, 1])



        


        #note we haven't implemented pf for r, seems pointless.
        #may be a better way to do this in the future as v little is actually changing in each loop, especiallly since θ and ζ will have a normal grid.
        #ie perhaps this could just be done for each r or something. Not sure this will be the slowdown tho.
        #TODO -> change the inputs to a single structure, this is cooked 
        create_local_basis!(Φ, S, dSr, dSθ, dSζ, ddSrr, ddSrθ, ddSrζ, ddSθθ, ddSθζ, ddSζζ, grids.θ.pf, grids.ζ.pf, dr, dθ, dζ)
        #create_local_basis!(Φ, S, dSr, dSθ, ddSrr, ddSrθ, ddSθθ, grids.θ.pf, n1, dr, dθ)


        #negatives for conjugate, will assume the phase factor is conjugate as well.
        create_local_basis!(Ψ, S, dSr, dSθ, dSζ, ddSrr, ddSrθ, ddSrζ, ddSθθ, ddSθζ, ddSζζ, -grids.θ.pf, -grids.ζ.pf, dr, dθ, dζ)



        for trialr in 1:4, trialθ in 1:4, trialζ in 1:4

            #display("got to here ok!")
            #may need a θN or something!
            right_ind = grid_to_index(i, j, k, trialr, trialθ, trialζ, grids)

            for testr in 1:4, testθ in 1:4, testζ in 1:4
                #display("testsf")
                #display(testsf)

                
                left_ind = grid_to_index(i, j, k, testr, testθ, testζ, grids)

                #only check for boundaries if this is true
                #no other i's can possibly give boundaries
                
                if i==1 || i==grids.r.N-1


                    if left_ind == right_ind && left_ind in boundary_inds

                        #rows[arr_count] = left_ind
                        #cols[arr_count] = right_ind
                        #Wdata[arr_count] = 1.0 + 0.0im
                        #Idata[arr_count] = 1.0 + 0.0im
                        #this is happening 4 times...
                        push!(rows, left_ind)
                        push!(cols, right_ind)
                        push!(Wdata, 0.25 + 0.0im)
                        push!(Idata, 0.25 + 0.0im)
                        
                        #arr_count += 1

                        #bounds_count += 1
                    
                    #otherwise the boundaries are set to zero, which for sparse matrices
                    #is the same as leaving blank.
                    elseif left_ind in boundary_inds
                        continue
                    elseif right_ind in boundary_inds
                        continue
                    #otherwise a regular case for these indicies.
                    else
                        #display("Are we getting stuck here?")
                        #rows[arr_count] = left_ind
                        #cols[arr_count] = right_ind
                        push!(rows, left_ind)
                        push!(cols, right_ind)
                        

                        #TODO -> should be almost the same, 
                        Wsum = @views gauss_integrate(Ψ[testr, testθ, testζ, :, :, :, :], Φ[trialr, trialθ, trialζ, :, :, :, :], W, wgr, wgθ, wgζ, jac, grids.r.gp, grids.θ.gp, grids.ζ.gp)


                        #Isum = @views gauss_integrate(Ψ[:, testr, :, testθ, :, testζ, :], Φ[:, trialr, :, trialθ, :, trialζ, :], I, wgr, wgθ, wgζ, jac, grids.r.gp, grids.θ.gp, grids.ζ.gp)
                        Isum = @views gauss_integrate(Ψ[testr, testθ, testζ, :, :, :, :], Φ[trialr, trialθ, trialζ, :, :, :, :], I, wgr, wgθ, wgζ, jac, grids.r.gp, grids.θ.gp, grids.ζ.gp)
                        #Wsum = 1
                        #Isum = 1
                        #display("perhaps?")

                        push!(Wdata, Wsum)
                        push!(Idata, Isum)
                    end
                else
                    
                    #rows[arr_count] = left_ind
                    #cols[arr_count] = right_ind
                    push!(rows, left_ind)
                    push!(cols, right_ind)
                        

                    #Wsum = @views gauss_integrate(Ψ[:, testr, :, testθ, :, testζ, :], Φ[:, trialr, :, trialθ, :, trialζ, :], W, wgr, wgθ, wgζ, jac, grids.r.gp, grids.θ.gp, grids.ζ.gp)
                    Wsum = @views gauss_integrate(Ψ[testr, testθ, testζ, :, :, :, :], Φ[trialr, trialθ, trialζ, :, :, :, :], W, wgr, wgθ, wgζ, jac, grids.r.gp, grids.θ.gp, grids.ζ.gp)


                    Isum = @views gauss_integrate(Ψ[testr, testθ, testζ, :, :, :, :], Φ[trialr, trialθ, trialζ, :, :, :, :], I, wgr, wgθ, wgζ, jac, grids.r.gp, grids.θ.gp, grids.ζ.gp)

                    #Wsum = 1
                    #Isum = 1
                    push!(Wdata, Wsum)
                    push!(Idata, Isum)
                    
                end
            end
        end

    end


    #maybe more consisnt for this function to return the rows and data as per parallal case.
    Wmat = sparse(rows, cols, Wdata)
    Imat = sparse(rows, cols, Idata)


    return Wmat, Imat
end
