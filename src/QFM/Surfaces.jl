

#Distaster zone - STAY OUT!


#will almost certainly need more than this
#hate the order of this, probably makes sense if we understood the theory!
struct QFMSurfaceT
    ρ :: Float64 #stick with rho for now, but perhaps this should just be s. v unclear exactly what this does. #think it is just a surface label, which we pick to be the edge.
    scos :: Array{Float64, 2}
    tsin :: Array{Float64, 2}
    ssin :: Array{Float64, 2}
    tcos :: Array{Float64, 2}
end


#better version!
struct SurfaceITPT
    M :: Int64
    N :: Int64
    scos_itp :: Array{BSplineKit.SplineInterpolations.SplineInterpolation, 2}
    tsin_itp :: Array{BSplineKit.SplineInterpolations.SplineInterpolation, 2}
end

#will acti a bit like a main function
#not sure what the inptus will be or the form of the output.
#this functions is more garbage than anything
#we need to srsly optimize, just reusing memory allocations will make a big difference.
function construct_surfaces(plist, qlist, sguesslist, prob)

    #qlist and plist are assumed to be the same size.

    #p is like a poloidal mode number
    #q is like a toroidal mode number
    #we do not use (m, n) here because these are the rational surfaces
    #that we are searching for, i.e p/q

    #still have no idea what these numbers even mean!
    MM = 4
    M = 24
    N = 8
    #MM = 4
    #M = 1
    #N = 2
    #sguess = 0.4

    #no idea what the bounding surfaces are about, going to ignore for now.
    #also, we have not implemented the jacobian for the action.
    #that will probably increase the accuracy.

    met = MetT()
    B = BFieldT()
    #not sure if this is how this is done.
    surfaces = QFMSurfaceT[]

    for i in 1:1:length(plist)
        #obvs will need to pass the other stuff in here somehow.
        #scos, tsin, ssin, tcos = action(plist[i], qlist[i], prob)
        scos, tsin, ssin, tcos = action2(plist[i], qlist[i], prob, met, B, MM, M, N, sguesslist[i])

        ρ = scos[1, 1] #surface label.
        push!(surfaces, QFMSurfaceT(ρ, scos, tsin, ssin, tcos))
        @printf("Found %d of %d surfaces.\n", i, length(plist))
    end

    #now we add the bounding surfaces, assumed to be at 0, 1
    #subject to chaneg obvs.

    #convert to r from flux.
    #ρ1 = sqrt(2*0.1)
    #ρ2 = sqrt(2*0.9)

    #push!(surfaces, straighten_boundary(ρ1, MM, M, N, prob))
    #push!(surfaces, straighten_boundary(ρ2, MM, M, N, prob))

    #unsure if we actually want to do this, 
    #as I don't think we want the surface objects tbh, just want the interpolant build from them
    return surfaces


end

#this seems an odd function to use
#seems like just picking a rational surface near the boundary and running the og alg would be better?
function straighten_boundary(ρ, MM, M, N, prob, tol=1e-9, niter=10)

    #need to start passing these around properly.
    #M = 8
    #N = 8

    nfft_θ = MM * M
    nfft_ζ = MM * N

    r = ones(1, nfft_θ, nfft_ζ) .* ρ
    ζ = LinRange(0, 2π, nfft_ζ+1)[1:end-1]
    θ = LinRange(0, 2π, nfft_θ+1)[1:end-1]

    ζarr = zeros(1, nfft_θ, nfft_ζ)

    #this is probably the wrong shape, needs to be the same shape as fout
    #from irfft2d.
    #can just set the size equal if needed.
    θarr = zeros(nfft_θ, nfft_ζ)

    for i in 1:nfft_θ, j in 1:nfft_ζ
        ζarr[1, i, j] = ζ[j]
    end

    λcos = zeros(M+1, 2*N+1)
    λsin = zeros(M+1, 2*N+1)
    iota = 0
    mlist = collect(range(0, M))
    nlist = [collect(0:N);collect(-N:-1)] 

    Bs = zeros(nfft_θ, nfft_ζ)
    Bθ = zeros(nfft_θ, nfft_ζ)
    Bζ = zeros(nfft_θ, nfft_ζ)

    met = MetT()
    B = BFieldT()

    for i in 1:niter
        iota_old = iota
        λcos_old = λcos
        λsin_old = λsin

        λreal = irfft2D(λcos, λsin, nfft_θ, nfft_ζ)

        for j in 1:nfft_θ, k in 1:nfft_ζ
            θarr[j, k] = θ[j] + λreal[j, k]

            prob.compute_met(met, ρ, θarr[j, k], ζarr[1, j, k], prob.geo.R0)
            compute_B!(B.B, met.J, prob.q, prob.isl, prob.isl2, ρ, θarr[j, k], ζarr[1, j, k])
            #zero chance this works.
            Bs[j, k] = B.B[1]
            Bθ[j, k] = B.B[2]
            Bζ[j, k] = B.B[3]
        end

        Bθ_o_Bζ = @. Bθ / Bζ

        cn, sn = rfft2D(Bθ_o_Bζ, M, N)

        iota = cn[1, 1]

        #mnarr = zerso(nfft_θ, nfft_ζ)
        #for j in 1:nfft_θ, k in 1:nfft_ζ
        #    mnarr[j, k] = 



        @. λcos[1, 2:end] = sn[1, 2:end] / (-mlist[1] * iota + nlist[2:end])

        @. λsin[1, 2:end] = cn[1, 2:end] / (+mlist[1] * iota - nlist[2:end])

        for j in 2:M, k in 1:N
            #    mnarr[j, k] = 
            λcos[j, k] = sn[j, k] / (-mlist[j] * iota + nlist[k])
            λsin[j, k] = cn[j, k] / (+mlist[j] * iota - nlist[k])
        end
        λcos[1, 1] = 0.0
        λsin[1, 1] = 0.0

        erriota = abs(iota - iota_old)
        errcn = maximum(abs.(λcos - λcos_old))
        errsn = maximum(abs.(λsin - λsin_old))

        if maximum([erriota, errcn, errsn]) < tol
            #return QFMSurfaceT(ρ, )
            break
        end


    end

    scos = zeros(size(λcos))
    scos[1, 1] = ρ
    ssin = zeros(size(λcos))

    return QFMSurfaceT(ρ, scos, λsin, ssin, λcos)

end




#shorthand for evaluateing our entire Interpolation matrix at once
#we will eventially want this to have the output as an input for memory reasons etc.
function itp_mat(surf_itp, r; deriv=0)

    #still not sure why it has this shape tbh!
    dim1 = surf_itp.M + 1
    dim2 = 2*surf_itp.N + 1

    
    #unsure if we need more than this!
    scos = zeros(dim1, dim2)
    tsin = zeros(dim1, dim2)
    if deriv == 0

        for i in 1:dim1
            for j in 1:dim2
                scos[i, j] = surf_itp.scos_itp[i, j](r)
                tsin[i, j] = surf_itp.tsin_itp[i, j](r)
            end
        end
    else
        for i in 1:dim1
            for j in 1:dim2
                scos[i, j] = (Derivative(deriv) * surf_itp.scos_itp[i, j])(r)
                tsin[i, j] = (Derivative(deriv) * surf_itp.tsin_itp[i, j])(r)
            end
        end 
    end
    return scos, tsin
end


function create_surf_itp(surfs)
    #surfs is an array of the surface struct, which needs some serious fixing.
    #so we need some rhosurfs for this???
    #hard copied from python, this are found during qfm construction I guess.
    #note that python also has the 'edges' problemo for another day I think.

    #hint that our data storage method is not ideal
    scosn = zeros((length(surfs), size(surfs[1].scos)[1], size(surfs[1].scos)[2]))
    tsinn = zeros((length(surfs), size(surfs[1].scos)[1], size(surfs[1].scos)[2]))

    #so I think we need to sort this shit.
    #so this sorts the surfaces by the ρ value, so that interpolation can work properly.
    #display(length(surfs))
    ρsurfs_unsort = zeros(length(surfs))

    for (i, surf) in enumerate(surfs)
        #won't bother with the others yet.

        ρsurfs_unsort[i] = surf.ρ
    end

    #may need to double check this is actually working properly!
    #looks to be sorting properly.
    perm = sortperm(ρsurfs_unsort)
    ρsurfs = ρsurfs_unsort[perm]
    surfs = surfs[perm]
    #so ρsurfs has to be in order
    #we can probably just have our inputs be in order, but 
    #we should also have some sorting here.
    #shouldn't be too hard.
    for (i, surf) in enumerate(surfs)
        #won't bother with the others yet.
        scosn[i, :, :] = surf.scos
        tsinn[i, :, :] = surf.tsin
        #ρsurfs[i] = surf.ρ
    end

    #ρsurfs = [0.6, 0.6111111111, 0.61538462, 0.61904762, 0.625]

    #these being hardcoded is really good.
    pqMpol = 24
    pqNtor = 8
    dim1 = pqMpol + 1
    dim2 = 2*pqNtor + 1

    #do I not need these here???
    #should check lol.
    mlist = collect(range(0, pqMpol))

    collect(-pqNtor:0)
    nlist = [collect(0:pqNtor);collect(-pqNtor:-1)] 
    #this is going to take a serious re-organise of the data.
    #we need to construct the grids to match the shape of scosn
    #so we need to have rho which matches the surfaces already.

    #this is probably what we want, but this requires that mlist and nlist are in order. that will require that we shuffle our data around. This will be v annoying.
    #scos_itp = interpolate((ρsurfs, mlist, nlist), scosn, (Gridded(Linear()), Gridded(Linear()), Gridded(Linear())))
    #not sure if this will do what we want. The surface labels have no correspondance to rho....
    #this is like obvs not going to work.

    #display(ρsurfs)

    #ok so this works, and outputs the 25x17 array we want
    #something like mpol+1 * 2npol or whatever.
    #not as robust, as it cannot handle stuff outside the edges.
    #so we either need to make this extrapolate
    #or haev some bounding box type stuff.
    #that is probably a non-coding issue. As the we will need to think
    #about how to deal with the non-chaotic region.
    #it may be best to only interpolate for specific r values.
    #or to find flux surfaces throughout the domain.
    #or do the boundaries which seems to be Zhisongs aproach.
    #ρsurfs = [ρsurfs[i] for i in 1:length(ρsurfs)]
    #scosn = [surfs[i].scos for i in 1:length(ρsurfs)]
    #tsinn = [surfs[i].tsin for i in 1:length(ρsurfs)]

    #display(ρsurfs)

    #this is satisfactory for now, but this function needs work
    #and we will need to compare more direclty with Zhisong
    #to make sure this form of interpolation is good enough!
    #scos_itp = interpolate((ρsurfs, ), scosn, Gridded(Linear()))
    #itp = interpolate(scosn, BSpline())
    #scos_itp = scale(itp, (ρsurfs, ))
    #scos_itp = scale(interpolate(scosn, BSpline(Linear())), ρsurfs )
    #tsin_itp = interpolate((ρsurfs, ), tsinn, Gridded(Linear()))

    #scos_itp = CubicSplineInterpolation(ρsurfs,  scosn)

    #lets try a matrix of interpolations. v silly but may work

    #itp_mat = Array(BSplineKit.SplineInterpolations, 25, 17)

    scos_itp = Array{BSplineKit.SplineInterpolations.SplineInterpolation}(undef, dim1, dim2)
    tsin_itp = Array{BSplineKit.SplineInterpolations.SplineInterpolation}(undef, dim1, dim2)

    #display(ρsurfs)
    #display(scosn[:,  10, 12])

    for i in 1:dim1
        for j in 1:dim2
            #display((i, j))
            #so interpolations cannot do unequal space other than linear.
            #fk me.
            #try the other package then.
            #it_mat[i, j] = cubic_spline_interpolation(ρsurfs, scosn[:, i, j])#, Gridded(Linear()))
            #we may not want to actually store this obkect
            #perhaps it would be better to store the approximate functions?
            #who knows tbh.
            #may want 3 to be a param, but for second deriv this is the minimum
            scos_itp[i, j] = interpolate(ρsurfs, scosn[:, i, j], BSplineOrder(3))
            tsin_itp[i, j] = interpolate(ρsurfs, tsinn[:, i, j], BSplineOrder(3))
        end
    end
    #tsin_itp = CubicSplineInterpolation(ρsurfs,  tsinn)

    #scos_itp = Spline1D(ρsurfs, scosn, k=1)


    #display(scos_itp)
    #display(tsin_itp(0.62))
    #ok, so to get this to work.
    #we want to interpolate θ, ζ (or their fourier values or whatever)
    #based on a single ρ value.
    #so we need to rewrite scosn to be an array of length ρ length
    #and each point contains the θ and ζ values.
    #we may need to flatten the scosn array in θ and ζ
    #scos_itp = interpolate((Base.OneTo(5), ρsurfs, ρsurfs), scosn, (NoInterp(), Gridded(Linear()), Gridded(Linear())))
    #, BSpline())
    #tsin_itp = interpolate(tsinn, BSpline())

    #display(scos_itp(0.62))

    #return SurfaceITPT(pqMpol, pqNtor, scos_itp, tsin_itp)
    return SurfaceITPT(pqMpol, pqNtor, scos_itp, tsin_itp)


end


#this function may not be needed tbh.
function convert_surf(surf)
    #change this name lol
    #takes the weird output form into a plotable form.

    #hard coded, yes very nice

    Nθ = 100
    θgrid = LinRange(0, 2*π, Nθ)
    ζ = 0.0

    #again, v stupid.
    scos = surf.scos
    ssin = surf.ssin
    tsin = surf.tsin
    tcos = surf.tcos

    α = zeros((Nθ, size(scos)[1], size(scos)[2]))

    pqMpol = 24
    pqNtor = 8
    mlist = collect(range(0, pqMpol))

    collect(-pqNtor:0)
    nlist = [collect(0:pqNtor);collect(-pqNtor:-1)] 

    for i in 1:Nθ
        for j in 1:size(scos)[1]
            for k in 1:size(scos)[2]
                α[i, j, k] = mlist[j] * θgrid[i] - nlist[k] * ζ
            end
        end
    end


    #α = @. 5 * θgrid - 8 * ζ



    cosα = cos.(α)
    sinα = sin.(α)

    s = zeros(Nθ)
    t = zeros(Nθ)

    for i in 1:Nθ

        for j in 1:size(scos)[1]
            for k in 1:size(scos)[2]
                s[i] += scos[j, k] * cosα[i, j, k]
                t[i] += θgrid[i] + tsin[j, k] * sinα[i, j, k]
            end
        end
    end

    return s, t

end


function farey_tree(N, q1, p1, q2, p2)

    #now N is a recursive depth.

    #think this is a stupid assertion.
    #won't help in weird cases,
    #also think it depends on the order.
    #@assert q2*p1 - q1*p2 == -1
    plist = [p1, p2]
    qlist = [q1, q2]

    farey_tree_recurs!(N, q1, p1, q2, p2, qlist, plist)

    #println(qlist)
    #println(plist)
    for i in eachindex(qlist)
        a = gcd(plist[i], qlist[i])
        if a != 1
            plist[i] ÷= a
            qlist[i] ÷= a
        end
    end

    #think we are going to assume that there are no double ups

    #unsure if this is automatically free from doubles??
    #now do some stuff with qlist and plist
    #so here want want to check there are not any doubles.
    #and reduce the gcd etc.
    #so this seems to work pretty well, however, this will only ever give surfaces between the islands, unsure how we will do them outside...
    return qlist, plist

end


function farey_tree_recurs!(N, q1, p1, q2, p2, qlist, plist)

    if N != 0
        N = N-1
        pnew = p1 + p2
        qnew = q1 + q2
        #worried about this tbh.
        push!(qlist, qnew)
        push!(plist, pnew)
        farey_tree_recurs!(N, qnew, pnew, q2, p2, qlist, plist)
        farey_tree_recurs!(N, q1, p1, qnew, pnew, qlist, plist)
        #return farey_tree_recurs(N, )
    end


end


#probably should just delete this.
function farey_tree_old(N, q1, p1, q2, p2)
    #N is the number of rationals we want to find,
    #eventually we will specify clustering etc.
    #first check the farey neighbouring condition thing

    @assert q1*p2 - q2*p1 == -1


    found = 0

    qlist = []
    plist = []

    iota1 = p1 / q1
    iota2 = p2 / q2

    #define the `parent` values.
    par_q1 = q1
    par_p1 = p1
    par_q2 = q2
    par_p2 = p2

    while found < N

        new_q = par_q1 + par_q2
        new_p = par_p1 + par_p2

        a = gcd(new_p, new_q)

        if a != 1
            new_p = new_p ÷ a
            new_q = new_q ÷ a
        end

        #check to make sure it is not already in the list.

        qinds = findall(x->x==new_q, qlist)
        pinds = findall(x->x==new_q, plist)

        if (any(qinds == pinds)) 

        else
            push!(plist, new_p)
            push!(qlist, new_q)
            found += 1
        end

        #i guess we will try to keep them central?
        #this is going to be real tricky
        #ideally this will be done recursively.
        new_iota = new_p / new_q



    end

    return qlist, plist


end
