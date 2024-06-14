
"""
Module for computing the island continuum.

"""
module IslandContinuum
#TODO

using MID.Misc
using MID.Geometry#; flux_toroidal_metric!

#think Elliptic is the go,
#maybe we should do the Elliptic; K, am, etc for just the funcs we need?
using Elliptic #may need to modify which of these we actually use.
using FFTW
using LinearAlgebra

#so looks like we used a few of the ol elliptic packages before, not ideal.

#so problem that needs to be understood, EllipticFunctions has Elliptk which does accept complex values
#but this is only needed because we are passing in a κ > 1, which is outside the domain, seems wrong
#previous code used Both EllipticFunctions and Elliptic, which seems like a bad idea.
#maybe some asserts on the input to stop that trapped and passing cases from cooking it?
#or a main function that calls both on correct parts of the list???


include("ContGeo.jl")
include("ContIsland.jl")


export trapped_continuum
export passing_continuum
export island_continuum

export ContIslandT

#so a few tricky parts to this.
#we probably want to do it in terms of r not ψ. Not clear how hard/possible that will be.
#also need to actually undertsand this.
#probably 
#I think end goal with this is to hopefully find the cut off for island amplitude A where this will
#cause damping, hopefully that will match our actual damping results.

#think we will just get it working first, then see if we can amke it work with our other functions etc and ideally swap to r.
#or should we swap everything to ψ... both are awful I think!

#I guess this is the time when we figure out how this actually works.


#maybe we just have a single function with sign=-1, 0, 1 so sign=0 means we are trapped?
#I think there is a lot of room to make this more similar to our other case,
#that will probably help it make sense.
#eg change the way I and W are computed, element by element not all at once
#change the way the fourier integration is done
#should probably just pass in grids tbh.
function passing_continuum(χlist, pmd, tmd, geo, isl, sign)

    #display(spectral_grid(pmd))
    nm, mlist, θ̄grid = spectral_grid(pmd) #bars are still a bit unclear tbh!
    nn, nlist, ζgrid = spectral_grid(tmd)

    ᾱ = collect(θ̄grid) .- collect(ζgrid)' ./ isl.q0 #not clear wot is going on here tbh.

    mlistbr = vec(repeat(mlist, 1, tmd.count)')
    nlistbr = repeat(nlist, pmd.count)


    #these just contain the index's for m and n, for simplified integration!
    midlist = vec(repeat(0:1:pmd.count-1, 1, tmd.count)')
    nidlist = repeat(0:1:tmd.count-1, pmd.count)
    #wont be doin the broadcasted lists! v confusing way to do fourier integrals.
    #but will start with them...
    #not exactly clear about the shapes and usage of this, but is doing the integration!
    diffm = midlist .- midlist'
    diffn = nidlist .- nidlist'

    #replicates pythons negative index behaviour
    #add 1 to match julia indexing, golly gosh this is a terrible solution.
    #maybe the way this is done in cflae would be clearer.
    diffm = mod.(diffm .+ nm, nm) .+ 1
    diffn = mod.(diffn .+ nn, nn) .+ 1

    ω2list = zeros(length(χlist), pmd.count * tmd.count)

    Imat = zeros(ComplexF64, size(diffm)...)
    Wmat = zeros(ComplexF64, size(diffm)...)

    ∇ψ̄2 = Matrix{Float64}(undef, nm, nn)
    #J = Matrix{Float64}(undef, nm, nn)
    #B2 = Matrix{Float64}(undef, nm, nn)
    #display(nm)
    B = B_field(Array{Float64}(undef, 3, nm, nn), Array{Float64}(undef, nm, nn))

    #shouldn't call this g!
    g = IslContMetT(zeros(3, 3, nm, nn), zeros(3, 3, nm, nn), Array{Float64}(undef, nm, nn))

    for (l, χ) in enumerate(χlist)
        #the way these functions are done is stupid af, may not be worth the effor to fix though!
        #this is comparable in time to python.
        #∇ψ̄2, J, B2 = ∇ψ̄2_J_B2_p(χ, ᾱ, geo, isl, metric, θgrid, ζgrid, sign)
        #∇ψ̄2_J_B2_p!(∇ψ̄2, J, B2, χ, ᾱ, geo, isl, metric, θgrid, ζgrid, sign)
        #display(isl.A)
        α = α_p(χ, ᾱ, isl)
        #display(α)
        ψ = ψ_p(χ, α, isl, sign)

        #shape of this is almost certain to not work, for slab case it won't matter
        #as that metric doesn't use theta
        θ = α .+ repeat(ζgrid, 1, nm)' ./isl.q0

        

        #display(α)
        #display(ζgrid ./isl.q0)
        #display(repeat(ζgrid, 1, nm)' ./isl.q0)
        #display(θ)
        
        #this will allow us to change which metric we use, but this is currently an ugly way of doing this!
        isl_cont_metric(g, ψ, θ, geo)
        compute_Bfield!(B, ψ, α, isl, g)

        #display(θ)

        ∇ψ̄2_p!(∇ψ̄2, χ, ψ, α, g, isl, sign)

        #display(∇ψ̄2)

        #seems like we need to scale this by J to match Zhisong, doesn't really make sense!
        #but we just dividing I, W by J does seem to effect the result, not clear why
        #not doing the division by nm and nn doesn't seem to effect things!
        #we may need to compare to some analytical solutions before we judge this!
        I = ∇ψ̄2 .* g.J ./ B.B2
        #display(I)
        
        #this is the part of W that does not include the derivative part
        #ie the extra factors of k,l,m,n
        #stupid name for this.
        W_no_par = I ./ g.J .^2 #Zhisong seems to have an extra J going around here, not sure why!
        #could be a mistake or we could have missed something!

        #return I

        q = 1/(ω_p(χ, isl, sign) + 1/isl.q0) # no idea where this comes from, we do compute q in Bfield, but differently, so this could be repeated?

        #also not sure what the hek this is
        mqn = mlistbr/q + nlistbr

        mqnmat = mqn .* mqn'

        #this does the same thing as python! not sure if this is the best way forward though
        #may want to consider plan_fft for speedup!

        Ifft = fft(I)
        Wfft = fft(W_no_par)
        
        #think we should know what this size is going to be ahead of time. ie m*n x m*n
        #looks like this is working as intended atm, not a very elegant solution though.
        #not sure this is actually doing what we want!
        #so this combined with how we define diffm and diffn are actually the problem.. :(
        #probably mostly because we don't know what this is supposed to be doing!
        #this might still be wrong...
        for i in 1:1:size(diffm)[1]
            for j in 1:1:size(diffn)[1]
                #println(diffm[i, i])
                #println(diffn[j, j])
                Imat[i, j] = Ifft[diffm[i, j], diffn[i, j]]
                Wmat[i, j] = Wfft[diffm[i, j], diffn[i, j]] * mqnmat[i, j] #need the mqn scaling still
            end
        end

        #need to factorise this or someshit first...
        #so big difference between here and python is that eigh in python assumes Hermitian
        #here we have to tell it, now the solving is similar.
        sol = eigen(Hermitian(Wmat), Hermitian(Imat))
        #print(sol.values)
        #need to normalise this probably.
        ω2list[l, :] = real.(sol.values) #before we can store these properly!


    end

    return ω2list

end


function island_W_and_I!(W, I, χ, θ̄grid, ζgrid, met, B, isl, geo, sign)

    #think we want to try and do this point by point.

    for (i, θ̄) in enumerate(θ̄grid), (j, ζ) in enumerate(ζgrid)

        #first we determine what ᾱ is.
        if sign == 0
            ᾱ = θ̄
        else
            ᾱ = θ̄ - ζ/isl.q0
        end
        #display(χ)
        #now we determine our original coordinate values, for computing the metric and B.
        α = compute_α(χ, ᾱ, isl, sign)

        
        ψ = compute_ψ(χ, α, ᾱ, isl, sign)

        #theta is actually different for passing vs trapped,
        #but zetagrid differences handle that.
        θ = α + ζ/isl.q0
        #display(θ)
        #could be the wrong ζ but met is indep of it so doesn't matter.
        flux_toroidal_metric!(met, ψ, θ, ζ, geo.R0)

        compute_island_B!(B, met, isl, ψ, α)
        #display(θ)
        #we want to compute 
        ∇χ2 = compute_∇χ2(ψ, α, met, isl)

        #display(∇χ2)

        #not clear wot mqn is, Think it is the parrallel gradient, will need to be handled
        #in the main section I think!
        #B will be annoying. Could reuse but that would be wasteful and we would have to pass in the q-profile
        #maybe better to just write a small B-field function specifically for the island.
        W[i, j] = ∇χ2 / (met.J * B.B2)

        I[i, j] = ∇χ2 * met.J / B.B2
    end
        

end


"""
So inputs assume a form of Φ = ∑ Φ_{m,n} exp(im(ANGLE1) - in(ANGLE2))

In both cases ANGLE2 is ζ, however in the trapped case, the domain of ζ has to be scaled by m0
this reflects that there are m0 islands ine ach tortoidal cross section.

ANGLE1 is more complicated. Outside the island, ANGLE1 is θ̄, but for inside the island 
ANGLE1 is ᾱ. It is not clear why this is done. It could be because ooutside θ̄ has a nicer interpretation, as it is almost a poloidal angle. While inside, ᾱ is an angle relative to the island centre.
#I think that is the answer, α is better inside, θ is better outside.

α = θ - ζ/q0 and the same is true for bars.

"""
#maybe region would be better than sign?
#sign = -1, 0, 1, 0 for trapped.
#so this looks like it is working now!
function island_continuum(χlist, pmd, tmd, geo, isl, sign)

    #in all cases inouts are taken to be θ̄ and ζ, in trapped case, θ̄ is equivalent to ᾱ.
    nθ, mlist, θ̄grid = spectral_grid(pmd) 
    nζ, nlist, ζgrid = spectral_grid(tmd)

    ᾱ = zeros(Float64, length(nθ), length(nζ))

    #trapped case
    if sign==0
        #inside the island the zeta grid is cooked, has to be extended
        #this accounts for there being m0 islands in a toroidal cross section
        #ie ζ is not single valued.
        ζgrid = LinRange(0, 2 * π * isl.m0 / tmd.incr, nζ+1)[1:end-1]
        #for the trapped case ᾱ is our poloidal angle, here we we expand it to the correct shape
        for i in 1:length(nθ), j in 1:length(nζ)
            ᾱ[i, j] = θ̄grid[i]
        end
    else
        for i in 1:length(nθ), j in 1:length(nζ)
            ᾱ[i, j] = θ̄grid[i] - ζgrid[j]/isl.q0
        end
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
    B = B_field(Array{Float64}(undef, 3, nθ, nζ), Array{Float64}(undef, nθ, nζ))

    B = IslandBT(zeros(3), 0.0)

    #metric elements need to be zeros, because they don't all get modified.
    #needs to change.
    met = MetT(zeros(3, 3), zeros(3, 3), zeros(3, 3, 3), zeros(3, 3, 3), 0.0, zeros(3))

    #seems odd to do this one but not the rest.
    #∇ψ̄2 = Matrix{Float64}(undef, nθ, nζ)

    p = plan_fft!(W, [1, 2])

    for (i, χ) in enumerate(χlist)

        island_W_and_I!(W, I, χ, θ̄grid, ζgrid, met, B, isl, geo, sign)

        #display(met.J)
        #display(B.B2)


        if sign == 0
            q = 1/(ω_t(χ, isl))
        else
            q = 1/(ω_p(χ, isl, sign) + 1/isl.q0)
        end
        #return I
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
                if sign == 0
                    gradfact = (m1/q + n1/isl.m0) *  (m2/q + n2/isl.m0)
                else 
                    gradfact = (m1/q + n1) *  (m2/q + n2)
                end
                #probably want a function for this!
                #last factor is a guess.
                #seems like one of those factors should be negative? for conjuagte?
                #but maybe that cancels with i^2?
                Wmat[left_ind, right_ind] = W[mind, nind] * gradfact
                Imat[left_ind, right_ind] = I[mind, nind] 
            end
        end
        #display(Imat)
        #sol = eig(Hermitian(Wmat), Hermitian(Imat))
        #ω2list[i, :] = real.(sol.values) 
                        
        ω2list[i, :] = real.(eigvals(Hermitian(Wmat), Hermitian(Imat)))
    end

    return ω2list

end


function trapped_continuum(χlist, pmd, tmd, geo, isl)

    #should be nθ and nζ to match everything else.
    nm, mlist, θ̄grid = spectral_grid(pmd) #bars are still a bit unclear tbh!
    nn, nlist, _ = spectral_grid(tmd)

    #inside the island the zeta grid is cooked, has to be extended
    ζgrid = LinRange(0, 2 * π * isl.m0 / tmd.incr, nn+1)[1:end-1]
    #display(θ̄grid)
    #display(ζgrid)
    #seems like we don't need the extra dimension for the trapped case???
    #this still does not make any sense. The nature of the grids is very unclear.
    #Really have no idea why the grids are the way they are lol.
    #also don't understand why we choose to use different angles for inside and outside.
    #not sure if this should be bar or not for θ, not sure if we 
    #need to undertsand this code fully tbh.
    ᾱ = θ̄grid .- ζgrid'/isl.q0 * 0

    display(size(ᾱ))


    #looks to be working as intended, is kinda slow though!
    #need to recompute nn and nm here I think. scope is awful! file is awful!
    mlistbr = vec(repeat(mlist, 1, tmd.count)')
    nlistbr = repeat(nlist, pmd.count)

    #these just contain the index's for m and n, for simplified integration!
    midlist = vec(repeat(0:1:pmd.count-1, 1, tmd.count)')
    nidlist = repeat(0:1:tmd.count-1, pmd.count)


    #not exactly clear about the shapes and usage of this, but is doing the integration!
    diffm = midlist .- midlist'
    diffn = nidlist .- nidlist'

    #display(size(diffm))
    #replicates pythons negative index behaviour
    #add 1 to match julia indexing, golly gosh this is a terrible solution.
    #this was the problem all along....
    diffm = mod.(diffm .+ nm, nm) .+ 1
    diffn = mod.(diffn .+ nn, nn) .+ 1


    ω2list = zeros(length(χlist), pmd.count * tmd.count) #not sure how big this bad boy needs to be, not sure how many evals are returned on each pass

    #I = Matrix{Float64}(undef, nm, nn)
    #W_no_par = Matrix{Float64}(undef, nm, nn)
    #Imat = Matrix{Float64}(undef, nm*nn, nm*nn)
    #Wmat = Matrix{Float64}(undef, nm*nn, nm*nn)
    Imat = zeros(ComplexF64, size(diffm)...)
    Wmat = zeros(ComplexF64, size(diffm)...)


    ∇ψ̄2 = Matrix{Float64}(undef, nm, nn)
    #J = Matrix{Float64}(undef, nm, nn)
    #B2 = Matrix{Float64}(undef, nm, nn)
    #display(nm)
    B = B_field(Array{Float64}(undef, 3, nm, nn), Array{Float64}(undef, nm, nn))

    #g = metric(Array{Float64}(undef, 3, 3, nm, nn), Array{Float64}(undef, 3, 3, nm, nn), Array{Float64}(undef, nm, nn))

    #metric elements need to be zeros, because they don't all get modified.
    g = IslContMetT(zeros(3, 3, nm, nn), zeros(3, 3, nm, nn), Array{Float64}(undef, nm, nn))

    #display(size(Imat))

    for (l, χ) in enumerate(χlist)
        #the way these functions are done is stupid af, may not be worth the effor to fix though!
        #println("Variable time")
        #this seems to be the biggest slow down compared to python.
        #@time ∇ψ̄2_J_B2_t!(∇ψ̄2, J, B2, χ, ᾱ, geo, isl, metric, θgrid, ζgrid)
        #display(χ)
        α = α_t(χ, ᾱ, isl)
        
        ψ = ψ_t(χ, α, ᾱ, isl)
        #display(ψ)
        #display(α)
        #got to be a better way to make this the shape we want!
        θ = α .+ repeat(ζgrid, 1, nm)' ./isl.q0
        #display(θ)
        #hate this lol
        isl_cont_metric(g, ψ, θ, geo)
        compute_Bfield!(B, ψ, α, isl, g)

        #display(g.gu)

        ∇ψ̄2_t!(∇ψ̄2, χ, ψ, α, g, isl)
        #display(∇ψ̄2)

        #display(∇ψ̄2)
        #testing if we actually need to multiply by j
        #multiply by J here is needed to match Zhisong, not certain it is correct!
        I = ∇ψ̄2 .* g.J ./ B.B2

        #display(g.J)
        #display(B.B2)

        #display(size(I))

        #this is the part of W that does not include the derivative part
        #ie the extra factors of k,l,m,n
        #stupid name for this.
        W_no_par = I ./ g.J .^2

        #display(W_no_par)

        q = 1/(ω_t(χ, isl)) # no idea where this comes from, we do compute q in Bfield, but differently, so this could be repeated?

        mqn = mlistbr/q + nlistbr/isl.m0

        mqnmat = mqn .* mqn'

        #return I
        #println("FTT time")
        #display(I)
        Ifft = fft(I, (1, 2))
        Wfft = fft(W_no_par)

        #display(Ifft)
       
        #display(Ifft)
        #think we should know what this size is going to be ahead of time. ie m*n x m*n
        #looks like this is working as intended atm, not a very elegant solution though.
        #not sure this is actually doing what we want!
        #so this combined with how we define diffm and diffn are actually the problem.. :(
        #probably mostly because we don't know what this is supposed to be doing!
        #println("loop time")
        #this is an absurdly expensive loop!
        #may need to try doing this withot diffm and diffn to see how it works.
        
        for i in 1:1:size(diffm)[1]
            for j in 1:1:size(diffn)[1]
                #println(diffm[i, i])
                #println(diffn[j, j])
                Imat[i, j] = Ifft[diffm[i, j], diffn[i, j]]
                Wmat[i, j] = Wfft[diffm[i, j], diffn[i, j]] * mqnmat[j, i] #need the mqn scaling still
            end
        end
        #display(diffm)
        #display(diffn)
        
        #alternative approach, may want to do this so it is more similar to our regualr code
        """
        for (i1, m1) in enumerate(mlist)
            for (i2, m2) in enumerate(mlist)
                for (j1, n1) in enumerate(nlist)
                    for (j2, n2) in enumerate(nlist)
                        #so this doesn't work, just appears to when ncount = 1
                        #so we need to get this to work!
                        mind = mod(m2-m1 + nm, nm) + 1
                        nind = mod(n2-n1 + nn, nn) + 1
                        #display(mind)
                        #display(nind)
                        Imat[j1 + (i1-1)*ncount, j2 + (i2-1)*ncount] += Ifft[mind, nind]
                        #not sure how to do the mqn matrix, can;t be bothered
                        #Wmat[i1 + (j1-1)*nm, i2 + (j2-1)*nm] = Wfft[mind, nind] * 
                    end
                end
            end
        end
        """
        
        #display(Imat)
        #println("Solve time")
        sol = eigen(Hermitian(Wmat), Hermitian(Imat))
        ω2list[l, :] = real.(sol.values) 

    end
    return ω2list



end

end