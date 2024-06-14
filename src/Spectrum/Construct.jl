"""
Constructs the two matrices using the WeakForm of the SAW governing equation. Uses Finite elements with cubic Hermite polynomials in r, and the fourier spectral method in θ and ζ. Returns two sparse matrices.

# Args
prob::ProblemT - Struct containing the functions and parameters that define the problem we are solving
grids::GridT - Grids to solve over.
"""
function construct(; prob::ProblemT, grids::GridsT)

    #island being nothing is stupid, either need to separate cases with and without island, or construct an empty (A=0) island if it is nothing.

    nθ, mlist, θgrid = spectral_grid(grids.pmd)
    nζ, nlist, ζgrid = spectral_grid(grids.tmd)


    #initialise the two structs.
    met = MetT(zeros(3, 3), zeros(3, 3), zeros(3, 3, 3), zeros(3, 3, 3), 0.0, zeros(3))
    B = BFieldT(zeros(3), zeros(3), zeros(3, 3), zeros(3, 3), 0.0, zeros(3))

    ξ, wg = gausslegendre(grids.rd.gp) #same as python!

    #gets the basis 
    H, dH, ddH = hermite_basis(ξ)

    #the trial function
    #4 is the number of Hermite shape functions
    #10 is Φ and all its relevant derivatives.
    #under new method we don't use the zeroth derivative, so these could be replaced with 9 
    Φ = zeros(ComplexF64, 9, 4, grids.rd.gp)
    #the test function.
    Ψ = zeros(ComplexF64, 9, 4, grids.rd.gp)   


    #generalised eval problem WΦ = ω^2 I Φ
    #Imat = zeros(ComplexF64, 2 * rd.grid_size * pmd.count * tmd.count, 2 * rd.grid_size * pmd.count * tmd.count)
    #Wmat = zeros(ComplexF64, 2 * rd.grid_size * pmd.count * tmd.count, 2 * rd.grid_size * pmd.count * tmd.count)

    #probably possible to know the size of this first.
    #can probably determine the maximum sized int we need based on size of matrix.
    #need to determine the size of these bad bois, they are eating up a bit of time.
    #this seems to have made a minimal difference, but is surely better practise right.
    arr_length = compute_length(grids.rd.N, grids.pmd.count, grids.tmd.count)

    arr_count = 1 #this may be wrong... gives error for v small matrix...

    rows = Array{Int64}(undef, arr_length) #ints
    cols = Array{Int64}(undef, arr_length) #ints
    Idata = Array{ComplexF64}(undef, arr_length)
    Wdata = Array{ComplexF64}(undef, arr_length)

    #rows = Array{Int64}(undef, 0) #ints
    #cols = Array{Int64}(undef, 0) #ints
    #Idata = Array{ComplexF64}(undef, 0)
    #Wdata = Array{ComplexF64}(undef, 0)


    #either need a condition in case m=0 or just an error message.
    boundary_inds = compute_boundary_inds(grids.rd.N, grids.pmd.count, grids.tmd.count, collect(mlist))
    #display(size(boundary_inds))
    #display(boundary_inds)

    I = zeros(ComplexF64, 9, 9, grids.rd.gp, nθ, nζ)
    W = zeros(ComplexF64, 9, 9, grids.rd.gp, nθ, nζ)


    #this gives a warning but seems to work perfectly
    #and incredibly efficiently!
    #this essentially defines a plane to fft I, (and W as they are the same size), which can be exectued
    #via p * I, this is done in place and seems to be mad efficient.
    p = plan_fft!(W, [4, 5])

    #Isum = 0.0 + 0.0im
    #Wsum = 0.0 + 0.0im

    bounds_count = 0

    #now we loop through the grid

    rgrid = construct_rgrid(grids)

    for i in 1:grids.rd.N-1

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
                mind = mod(k1-k2 + nθ, nθ) + 1
                nind = mod(l1-l2 + nζ, nζ) + 1


                for trialsf in 1:4

                    right_ind = grid_to_index(i, k1, l1, trialsf, grids.pmd.count, grids.tmd.count)

                    for testsf in 1:4
                        #display("testsf")
                        #display(testsf)

                        
                        left_ind = grid_to_index(i, k2, l2, testsf, grids.pmd.count, grids.tmd.count)

                        #only check for boundaries if this is true
                        #no other i's can possibly give boundaries
                        if i==1 || i==grids.rd.N-1


                            if left_ind == right_ind && left_ind in boundary_inds

                                rows[arr_count] = left_ind
                                cols[arr_count] = right_ind
                                Wdata[arr_count] = 1.0 + 0.0im
                                Idata[arr_count] = 1.0 + 0.0im
                                
                                arr_count += 1

                                #bounds_count += 1
                            
                            #otherwise the boundaries are set to zero, which for sparse matrices
                            #is the same as leaving blank.
                            elseif left_ind in boundary_inds
                                continue
                            elseif right_ind in boundary_inds
                                continue
                            #otherwise a regular case for these indicies.
                            else
                                rows[arr_count] = left_ind
                                cols[arr_count] = right_ind
                                

                                Wsum = @views gauss_integrate(Ψ[:, testsf, :], Φ[:, trialsf, :], W[:, :, :, mind, nind], wg, jac, grids.rd.gp)


                                Isum = @views gauss_integrate(Ψ[:, testsf, :], Φ[:, trialsf, :], I[:, :, :, mind, nind], wg, jac, grids.rd.gp)

                                Wdata[arr_count] = Wsum
                                Idata[arr_count] = Isum
                                
                                arr_count += 1
                            end
                        else
                            
                            rows[arr_count] = left_ind
                            cols[arr_count] = right_ind
                                

                            Wsum = @views gauss_integrate(Ψ[:, testsf, :], Φ[:, trialsf, :], W[:, :, :, mind, nind], wg, jac, grids.rd.gp)


                            Isum = @views gauss_integrate(Ψ[:, testsf, :], Φ[:, trialsf, :], I[:, :, :, mind, nind], wg, jac, grids.rd.gp)

                            Wdata[arr_count] = Wsum
                            Idata[arr_count] = Isum
                            
                            arr_count += 1
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