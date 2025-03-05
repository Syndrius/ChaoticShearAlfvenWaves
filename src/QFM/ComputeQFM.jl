



#probs that will defs occur
#interpolation outside given region.
#need to add boundaries I guess
#interpolation is not currently working...


#think we will probably want non-continuum version eventually?
#no fkn idea what inputs this will need.
#ok so we can get through this without it all breaking!
#we are getting complete gibberish though!
#still, good progress!
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
    rgrid = LinRange(1.1095, 1.112, 10)

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

        W_and_I!(W, I, tor_met, tor_B, qfm_met, qfm_B, prob, [r], θgrid, ζgrid, tm, surf_itp, CT)

        #display(CT.JM)
        #display(tor_B.B)

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
        
        vals = eigvals(Hermitian(Wmat), Hermitian(Imat))


        ωlist[i, :] = prob.geo.R0 * sqrt.(abs.(vals))
    end

    #should probably return an evals Obj.
    return ωlist


end




#is W and I a stupid name?
#this kind of needs to be here because of the qfm strutcs,
#ideally this would not be the case.
#I guess we could store them in the structure module???
#seems sub-optimal.
#but this is not a good place to keep this function
function W_and_I!(W::Array{ComplexF64, 5}, I::Array{ComplexF64, 5}, tor_met::MetT, tor_B::BFieldT, qfm_met::MetT, qfm_B::BFieldT,   prob::ProblemT, r::Array{Float64}, θ::AbstractArray, ζ::AbstractArray, tm::TM, surf::SurfaceITPT, CT::CoordTsfmT)

    #compute the density.
    n = prob.dens.(r) :: Array{Float64}
    #TODO
    #ωcap2 = ω_cap2.(r) :: Array{Float64}
    ωcap2 = zeros(length(r))

    for k=1:1:length(ζ), j=1:1:length(θ), i=1:1:length(r)

        #so I am unsure if the new coords are ever actually used???
        coord_transform!(r[i], θ[j], ζ[k], CT, surf)
        #display(CT.JM)

        #compute the original metric
        prob.compute_met(tor_met, r[i], θ[j], ζ[k], prob.geo.R0)
        #and original B field.
        compute_B!(tor_B, tor_met, prob.q, prob.isl, prob.isl2, r[i], θ[j], ζ[k])

        #transform the metric
        met_transform!(tor_met, qfm_met, CT)

        #transform the B field
        B_transform!(tor_B, qfm_B, qfm_met, CT)


        #now we compute the weakform in the usual way.
        #computes the matrix D.
        WeakForm.compute_D!(qfm_met, qfm_B, tm.D)

        #compute the W matrix
        @views WeakForm.compute_W!(W[:, :, i, j, k], qfm_met, qfm_B, n[i], ωcap2[i], tm)

        #compute the I matrix
        @views WeakForm.compute_I!(I[:, :, i, j, k], qfm_met, qfm_B, n[i], prob.flr, tm.D, tm.F)

    end

end