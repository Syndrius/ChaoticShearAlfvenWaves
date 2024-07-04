

function continuum(prob::ProblemT, grids::ContGridsT, perN=true::Bool)

    #we can do this because usual cases (i.e. no island) have no toroidal coupling
    #so we can compute the continuum for each n individually
    if perN
        _, _, _, _, _, nlist, _ = instantiate_grids(grids)

        ωlist = zeros(grids.rN, grids.θ.count, grids.ζ.count)

        for (i, n) in enumerate(nlist)
            temp_ζgrid = init_sm_grid(start=n, count=1)
            temp_grids = init_grids(grids.rN, grids.θ, temp_ζgrid)

            ωlist[:, :, i] = compute_continuum(prob, temp_grids)
        end
    else
        #are we ever actually going to want this??
        ωlist = compute_continuum(prob, grids)
    end

    return ωlist

end


""" 
Finds the continuum by solving the second order derivatives on each flux surface individually. Much faster than reconstructing the continuum from the full solver, but cannot handle islands, resistivity and won't find any global modes.

# Args
prob::ProblemT - Struct containing the functions and parameters that define the problem we are solving
grids::GridT - Grids to solve over.
"""
function compute_continuum(prob::ProblemT, grids::ContGridsT)

    #don't think this will work with fem the whole time.

    #maybe should assert that the rgrid doesn't start from zero, as fem case avoids the zero,
    #cont case does not.

    #maybe easier to just assume that it will start from 0???
    #for contiinuum case we ignore the clustered region.
    #rlist = collect(LinRange(0, 1, grids.r.N))[2:end]

    rgrid, Nθ, mlist, θgrid, Nζ, nlist, ζgrid = instantiate_grids(grids)

    

    if prob.isl.A != 0.0
        display("Continuum can't handle islands you goose.")
        #exit()
    end

    if prob.δ != 0.0
        display("Continuum can't handle resisitivity you goose.")
        #exit()
    end

    #nθ, mlist, θgrid = spectral_grid(grids.pmd)
    #nζ, nlist, ζgrid = spectral_grid(grids.tmd)

    met = MetT(zeros(3, 3), zeros(3, 3), zeros(3, 3, 3), zeros(3, 3, 3), 0.0, zeros(3))
    B = BFieldT(zeros(3), zeros(3), zeros(3, 3), zeros(3, 3), 0.0, zeros(3))



    mat_size = matrix_dim(grids)

    #efuncs are useless in this case I think.
    #ϕlist = zeros(matrix_dim, length(rlist), pmd.count, tmd.count)
    #don't think we actually want to do this anymore!
    #probably would be better to convert this straight to the modes asap. 
    #are we ever going to not want that?
    ωlist = zeros(grids.rN, mat_size)

    I = zeros(ComplexF64, 9, 9, 1, Nθ, Nζ)
    W = zeros(ComplexF64, 9, 9, 1, Nθ, Nζ)

    #plan wont work the same as we are changing the size of I.
    #could probably still be done tbh! would be more complicated and maybe not worth it tbh!
    #p = plan_fft!(I, [4, 5])

    #now we loop through the grid

    for (i, r) in enumerate(rgrid) 


        #compute_met(met, [r], θgrid, ζgrid, R0)
        #compute_B(B, met, [r], θgrid, ζgrid, isl) 
        #current term is not needed for continuum
        #I, Tl, Tj, Idamp = I_and_W(B, met, dens, [r], θgrid, ζgrid, 0.0) 

        #I_and_W!(I, W, B, q_profile, met, compute_met, dens, [r], θgrid, ζgrid, 0.0, isl, R0)
        #does heaps of extra computation, but the continuum is quite efficeint so does not matter.
        W_and_I!(W, I, met, B, prob, [r], θgrid, ζgrid)

        #For the continuum we extract the relevant second derivative parts
        Icont = I[[1], [1], :, :, :]
        #this excludes ϕ_ss which should be part, but I think this only occurs when B^s≠0, which would break the flux surfaces.
        Wcont = W[[1, 5, 6], [1, 5, 6], :, :, :]


        Ifft = fft(Icont, (4, 5))
        Wfft = fft(Wcont, (4, 5))
        
        Imat = zeros(ComplexF64, mat_size, mat_size)
        Wmat = zeros(ComplexF64, mat_size, mat_size)


        #loop over the fourier components of the trial function
        for (k1,m1) in enumerate(mlist), (l1, n1) in enumerate(nlist)
            #left_ind = cont_grid_to_index(k1, l1, grids.ζ.count)
            left_ind = grid_to_index(k1, l1, grids)

            #equivalent to create local basis. Much simpler in contninuum case.
            Φ = [1, m1*1im, n1*1im]
            
            for (k2, m2) in enumerate(mlist), (l2, n2) in enumerate(nlist)
                #right_ind = cont_grid_to_index(k2, l2, grids.ζ.count)
                right_ind = grid_to_index(k2, l2, grids)
                
                Ψ = [1, -m2*1im, -n2*1im]

                #mind = mod(m1-m2 + nθ, nθ) + 1
                #nind = mod(n1-n2 + nζ, nζ) + 1
                mind = mod(k1-k2 + Nθ, Nθ) + 1
                nind = mod(l1-l2 + Nζ, Nζ) + 1

                
                

                #silly matrix dims, but keeps it consistent.
                Imat[left_ind, right_ind] += Ifft[1, 1, 1, mind, nind]

                #contract over the different derivatives.
                for j in 1:3
                    for k in 1:3
                        Wmat[left_ind, right_ind] += Ψ[j] * Wfft[j, k, 1, mind, nind] * Φ[k]
                    end
                         
                end
            end
        end

        #not sure if it makes sense to store the efuncs, could be worth looking at and comparing!
        #we will need to make another recosntruct function
        #think efuncs in this case are garbage
        vals = eigvals(Hermitian(Wmat), Hermitian(Imat))

        #ideally this would be in a better form I think!
        #imaginary part is basically zero.
        #I think the phi here is actually uselss as it doesn't contain the radial structure
        #a bit unclear though
        #ϕlist[:, i, :, :] = real.(cont_reconstruct_phi(funcs, matrix_dim, pmd.count, tmd.count))

        #abs only needed for cyl limit when very close to zero
        #might also be nice to somehow label each mode 
        #not sure how that will work, may need to use the efuncs.
        #would probbaly be best to use a try catch and just give a warning that some vals are negative?
        #or return the non-square root, so we can make sure its only very small numbers that are negative.
        #this reconstructs the mode form we want hopefully!
        #could probably do this after all are computed perhaps?
        #for j in 1:1:length(vals)

            #θind, ζind = index_to_grid(j, grids)
            #ωlist[i, θind, ζind] = prob.geo.R0 * sqrt.(abs.(vals[j]))
        #end
        ωlist[i, :] = prob.geo.R0 * sqrt.(abs.(vals))

    end

    
    return ωlist
end

