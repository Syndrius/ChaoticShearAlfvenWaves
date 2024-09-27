

"""
    island_continuum(χlist::Array{Float64}, pmd::Inputs.SMGridDataT, tmd::Inputs.SMGridDataT, geo::GeoParamsT, isl::IslandT, sign::Int64)

Computes the continuum with a magnetic island chain. Exact for case with linear rotational transform as outlined by Qu and Hole. Can also approximated the local continuum around the island via a linear taylor expansion of the rotational transform.


So inputs assume a form of Φ = ∑ Φ_{m,n} exp(im(ANGLE1) - in(ANGLE2))

In both cases ANGLE2 is ζ, however in the trapped case, the domain of ζ has to be scaled by m0
this reflects that there are m0 islands ine ach tortoidal cross section.

ANGLE1 is more complicated. Outside the island, ANGLE1 is θ̄, but for inside the island 
ANGLE1 is ᾱ. It is not clear why this is done. It could be because ooutside θ̄ has a nicer interpretation, as it is almost a poloidal angle. While inside, ᾱ is an angle relative to the island centre.
#I think that is the answer, α is better inside, θ is better outside.

α = θ - ζ/q0 and the same is true for bars.

# Args
χlist::Array{Float64} - List of energies χ to consider.
pmd::ModeDataT - poloidal mode data ie modes for ANGLE1.
tmd::ModeDataT - toroidal mode data ie modes for ANGLE2.
geo::GeoParamsT - struct storing the geometrical parameters, i.e major radius.
isl::IslandT - expanded struct storing the island paramaters needed for the continuum calculation.
sign::Int64 - Sign of the particles, ±1 for passing on either side, or 0 for trapped particles.
"""
function island_continuum(χlist::Array{Float64}, pmd::SMGridDataT, tmd::SMGridDataT, geo::GeoParamsT, isl::IslandT, sign::Int64)

    #in all cases inouts are taken to be θ̄ and ζ, in trapped case, θ̄ is equivalent to ᾱ.
    nθ, mlist, θ̄grid = Structures.sm_grid(pmd) 
    nζ, nlist, ζgrid = Structures.sm_grid(tmd)

    if sign == 0
        ζgrid = range(0, 2*π / tmd.incr * isl.m0, nζ+1)[1:end-1]
    end

    
    ω2list = zeros(length(χlist), pmd.count * tmd.count)

    #these are the matrices to be solved.
    #not sure if we need to reset them to zero each time, 
    #think not.
    Wmat = zeros(ComplexF64, pmd.count * tmd.count, pmd.count * tmd.count)
    Imat = zeros(ComplexF64, pmd.count * tmd.count, pmd.count * tmd.count)

    #should change to nθ and nζ tbh! would match other code.
    W = zeros(ComplexF64, nθ, nζ)
    I = zeros(ComplexF64, nθ, nζ)

    #needs to change

    B = BFieldT(zeros(3), zeros(3), zeros(3, 3), zeros(3, 3), 0.0, zeros(3))

    #metric elements need to be zeros, because they don't all get modified.
    #needs to change.
    met = MetT(zeros(3, 3), zeros(3, 3), zeros(3, 3, 3), zeros(3, 3, 3), 0.0, zeros(3))


    p = plan_fft!(W, [1, 2])

    for (i, χ) in enumerate(χlist)

        island_W_and_I!(W, I, χ, θ̄grid, ζgrid, met, B, isl, geo, sign)


        #effective q-profile is different for trapped or passing.
        #could do some + (abs(sign)+1) in front of final term if we want.
        if sign == 0
            q = 1/(compute_Ω(χ, isl, sign))
        else
            q = 1/(compute_Ω(χ, isl, sign) + 1/isl.q0)
        end

        #display(I)
        #fft W and I
        p * W
        p * I


        for (k1, m1) in enumerate(mlist), (l1, n1) in enumerate(nlist)
            left_ind = l1 + (k1-1) * tmd.count
            for (k2, m2) in enumerate(mlist), (l2, n2) in enumerate(nlist)

                right_ind = l2 + (k2-1) * tmd.count
                mind = mod(k1-k2 + nθ, nθ) + 1
                nind = mod(l1-l2 + nζ, nζ) + 1

                #probably a wildly inefficient way to do this!
                #scaling factor from the parallal gradient
                #different for trapped vs passing particles.
                if sign == 0
                    gradfact = (m1/q + n1/isl.m0) *  (m2/q + n2/isl.m0)
                else 
                    gradfact = (m1/q + n1) *  (m2/q + n2)
                end

                Wmat[left_ind, right_ind] = W[mind, nind] * gradfact
                Imat[left_ind, right_ind] = I[mind, nind] 
            end
        end

        #sol = eig(Hermitian(Wmat), Hermitian(Imat))
        #ω2list[i, :] = real.(sol.values) 
                        
        ω2list[i, :] = real.(eigvals(Hermitian(Wmat), Hermitian(Imat)))
    end

    return ω2list

end
