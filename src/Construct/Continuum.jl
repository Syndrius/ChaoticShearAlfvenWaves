
""" 
    continuum(prob::ProblemT, grids::ContGridsT)

Finds the continuum by solving the second order derivatives on each flux surface individually. Much faster than reconstructing the continuum from the full solver, but cannot handle islands, resistivity and won't find any global modes.
"""
function continuum(prob::ProblemT, grids::ContGridsT)

    #don't think this will work with fem the whole time.

    #maybe should assert that the rgrid doesn't start from zero, as fem case avoids the zero,
    #cont case does not.


    rgrid, θgrid, ζgrid = inst_grids(grids)

    Nθ = length(θgrid)
    Nζ = length(ζgrid)

    mlist = mode_list(grids.θ)
    nlist = mode_list(grids.ζ)

    

    if prob.isl.A != 0.0
        display("Continuum can't handle islands you goose.")
        #exit()
    end

    if prob.flr.δ != 0.0
        display("Continuum can't handle resisitivity you goose.")
        #exit()
    end

    #nθ, mlist, θgrid = spectral_grid(grids.pmd)
    #nζ, nlist, ζgrid = spectral_grid(grids.tmd)

    met = MetT()
    B = BFieldT()


    mat_size = matrix_size(grids)

    #efuncs are useless in this case I think.
    #ϕlist = zeros(matrix_dim, length(rlist), pmd.count, tmd.count)
    #don't think we actually want to do this anymore!
    #probably would be better to convert this straight to the modes asap. 
    #are we ever going to not want that?
    ωlist = zeros(grids.r.N, mat_size)

    #guess we should make a function for this like the other cases.
    I = zeros(ComplexF64, 9, 9, 1, Nθ, Nζ)
    W = zeros(ComplexF64, 9, 9, 1, Nθ, Nζ)

    #plan wont work the same as we are changing the size of I.
    #could probably still be done tbh! would be more complicated and maybe not worth it tbh!
    #p = plan_fft!(I, [4, 5])

    tm = TM()

    #now we loop through the grid

    for (i, r) in enumerate(rgrid) 


        #compute_met(met, [r], θgrid, ζgrid, R0)
        #compute_B(B, met, [r], θgrid, ζgrid, isl) 
        #current term is not needed for continuum
        #I, Tl, Tj, Idamp = I_and_W(B, met, dens, [r], θgrid, ζgrid, 0.0) 

        #I_and_W!(I, W, B, q_profile, met, compute_met, dens, [r], θgrid, ζgrid, 0.0, isl, R0)
        #does heaps of extra computation, but the continuum is quite efficeint so does not matter.
        W_and_I!(W, I, met, B, prob, [r], θgrid, ζgrid, tm)

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

    #ideally this would return an evals like structure instead of this
    
    return ωlist
end




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
