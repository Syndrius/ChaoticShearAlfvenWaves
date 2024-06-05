

"""
Constructs the two matrices and solves. Can either solve the full spectrum via inbuilt solving (slow), or use a shift and invert to find nev amount of eigenvalues closest to σ (fast).

# Args
prob::ProblemT - Struct containing the functions and parameters that define the problem we are solving
grids::GridT - Grids to solve over.
efuncs::Bool - Return eigenfunctions with values.
σ::Float64=0.0 - Find nev nearest evals to σ when solving with arpack.
reconstruct::Bool - Whether to reconstruct the eigenfunctions into 3d.
full_spectrum::Bool - Whether to solve for the full spectrum with inbuilt solver (slow) or use arpack (fast).
nev::Int64 - Number of eigenvalues to solve for if using Arpack.
"""
function analytical_construct_and_solve(; prob::ProblemT, grids::GridsT, efuncs=true::Bool, σ=0.0::Float64, reconstruct=true::Bool, full_spectrum=false::Bool, nev=20::Int64)

    W, I = analytical_construct(prob=prob, grids=grids)
    
    if full_spectrum 
        if prob.δ == 0.0
            ω, ϕ = full_spectrum_solve(Wmat=W, Imat=I, grids=grids, efuncs=efuncs, reconstruct=reconstruct, resistivity=false, R0=prob.geo.R0)
        else
            ω, ϕ = full_spectrum_solve(Wmat=W, Imat=I, grids=grids, efuncs=efuncs, reconstruct=reconstruct, resistivity=true, R0=prob.geo.R0)
        end
    else
        ω, ϕ = arpack_solve(Wmat=W, Imat=I, grids=grids, efuncs=efuncs, nev=nev, σ=σ, reconstruct=reconstruct, R0=prob.geo.R0)
    end
end




"""
Constructs the two matrices using the WeakForm of the SAW governing equation. Uses Finite elements with cubic Hermite polynomials in r, and the fourier spectral method in θ and ζ. Returns two sparse matrices.

# Args
prob::ProblemT - Struct containing the functions and parameters that define the problem we are solving
grids::GridT - Grids to solve over.
"""
function analytical_construct(; prob::ProblemT, grids::GridsT)
    
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
    Φ = zeros(ComplexF64, 10, 4, grids.rd.gp)
    #the test function.
    Ψ = zeros(ComplexF64, 10, 4, grids.rd.gp)   


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

    I = zeros(ComplexF64, 10, 10, grids.rd.gp, nθ, nζ)
    W = zeros(ComplexF64, 10, 10, grids.rd.gp, nθ, nζ)


    #this gives a warning but seems to work perfectly
    #and incredibly efficiently!
    #this essentially defines a plane to fft I, (and W as they are the same size), which can be exectued
    #via p * I, this is done in place and seems to be mad efficient.
    p = plan_fft!(I, [4, 5])

    #Isum = 0.0 + 0.0im
    #Wsum = 0.0 + 0.0im

    bounds_count = 0

    #now we loop through the grid

    rgrid = construct_rgrid(grids)

    for i in 1:grids.rd.N-1

        r, dr = local_to_global(i, ξ, rgrid)

        jac = dr/2 #following thesis!


        #I_and_W!(I, W, B, q_profile, met, compute_met, dens, r, θgrid, ζgrid, δ, isl, R0)

        #W_and_I!(W, I, met, B, prob, r, θgrid, ζgrid)
        analytical_W_and_I!(W, I, met, B, prob, r, θgrid, ζgrid)
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
        for (k1,m1) in enumerate(mlist)
            for (l1, n1) in enumerate(nlist)

                create_local_basis!(Φ, H, dH, ddH, m1, n1, jac)

                for (k2, m2) in enumerate(mlist)
                
                    for (l2, n2) in enumerate(nlist)

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
                                        

                                        Wsum = @views gauss_integrate_for_big(Ψ[:, testsf, :], Φ[:, trialsf, :], W[:, :, :, mind, nind], wg, jac, grids.rd.gp)


                                        Isum = @views gauss_integrate_for_big(Ψ[:, testsf, :], Φ[:, trialsf, :], I[:, :, :, mind, nind], wg, jac, grids.rd.gp)

                                        Wdata[arr_count] = Wsum
                                        Idata[arr_count] = Isum
                                        
                                        arr_count += 1
                                    end
                                else
                                    
                                    rows[arr_count] = left_ind
                                    cols[arr_count] = right_ind
                                        

                                    Wsum = @views gauss_integrate_for_big(Ψ[:, testsf, :], Φ[:, trialsf, :], W[:, :, :, mind, nind], wg, jac, grids.rd.gp)


                                    Isum = @views gauss_integrate_for_big(Ψ[:, testsf, :], Φ[:, trialsf, :], I[:, :, :, mind, nind], wg, jac, grids.rd.gp)

                                    Wdata[arr_count] = Wsum
                                    Idata[arr_count] = Isum
                                    
                                    arr_count += 1
                                end
                            end
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




function analytical_W_and_I!(W::Array{ComplexF64, 5}, I::Array{ComplexF64, 5}, met::MetT, B::BFieldT, prob::ProblemT, r:: Array{Float64}, θ::LinRange{Float64, Int64}, ζ::LinRange{Float64, Int64})

    #define toroidal versions.
    W_tor = zeros(ComplexF64, size(W))
    I_tor = zeros(ComplexF64, size(I))
    met_tor = MetT(zeros(3, 3), zeros(3, 3), zeros(3, 3, 3), zeros(3, 3, 3), 0.0, zeros(3))
    B_tor = BFieldT(zeros(3), zeros(3), zeros(3, 3), zeros(3, 3), 0.0, zeros(3))

    n = prob.dens.(r) :: Array{Float64}

    for k=1:1:length(ζ), j=1:1:length(θ), i=1:1:length(r)

        #compute the metric
        #if we do this, it matches almost perfectly with the results from Axels equation, assuming his equation is in the correct form,
        #However, using the form in his paper doesn't seem to give his results, this requires changing some signs,
        #this then gives us the results in his paper, but then our code does not match anymore.
        #this is all a bit fked.
        #also the no_delta part seems inconsistent with Axel's paper, but is clearly stated in the appendix of the equation paper, ie where they define the metric elements.
        #doing the normal calculation with test metric, we get a fair bit closer for R0=10
        #but still a bit off, then with R0=20 its almost spot on.
        cylindrical_metric!(met, r[i], θ[j], ζ[k], prob.geo.R0)
        #test_metric!(met_tor, r[i], θ[j], ζ[k], prob.geo.R0)
        toroidal_metric!(met_tor, r[i], θ[j], ζ[k], prob.geo.R0)
        #no_delta_metric!(met_tor, r[i], θ[j], ζ[k], prob.geo.R0)

        #compute the magnetic field.
        compute_B!(B, met, prob.q, prob.isl, r[i], θ[j], ζ[k])
        compute_B!(B_tor, met_tor, prob.q, prob.isl, r[i], θ[j], ζ[k])
        

        #Tj, the current term, is always computed under the cylindrical approxmiation
        #so we need to split up our computation of W.
        Tl = WeakForm.new_compute_Tl(met, B) .* met.J
        Tl_tor = WeakForm.new_compute_Tl(met_tor, B_tor) .* met_tor.J
        Tj = WeakForm.new_compute_Tj(met, B) .* met.J .* WeakForm.new_jparonB(met, B) ./ 2
        #@views new_compute_W!(W[:, :, i, j, k], met, B)
        

        #assign the total value of W.
        W[1:9, 1:9, i, j, k] .= Tl .- Tj

        #replace the double radial derivative terms with the toroidal part.
        W[1, 1, i, j, k] = Tl_tor[1, 1] - Tj[1, 1]
        W[4, 1:9, i, j, k] .= Tl_tor[4, :] .- Tj[4, :]
        W[1:9, 4, i, j, k] .= Tl_tor[:, 4] .- Tj[:, 4]

        #same for I.
        @views WeakForm.compute_I!(I[:, :, i, j, k], met, B, n[i], prob.δ)
        @views WeakForm.compute_I!(I_tor[:, :, i, j, k], met_tor, B_tor, n[i], prob.δ)
        I[1, 1, i, j, k] = I_tor[1, 1, i, j, k]
        I[4, :, i, j, k] .= I_tor[4, :, i, j, k]
        I[:, 4, i, j, k] .= I_tor[:, 4, i, j, k]



    end

end