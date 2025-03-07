

#Distaster zone - STAY OUT!


#will almost certainly need more than this
struct QFMSurfaceT
    ρ :: Float64 #stick with rho for now, but perhaps this should just be s. v unclear exactly what this does. #think it is just a surface label, which we pick to be the edge.
    scos
    tsin
    ssin
    tcos
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
function construct_surfaces(plist, qlist, prob)

    #qlist and plist are assumed to be the same size.

    #p is like a poloidal mode number
    #q is like a toroidal mode number
    #we do not use (m, n) here because these are the rational surfaces
    #that we are searching for, i.e p/q


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
        scos, tsin, ssin, tcos = action2(plist[i], qlist[i], prob, met, B)

        ρ = scos[1, 1] #surface label.
        push!(surfaces, QFMSurfaceT(ρ, scos, tsin, ssin, tcos))
    end

    #unsure if we actually want to do this, 
    #as I don't think we want the surface objects tbh, just want the interpolant build from them
    return surfaces




end

#this should probably be in plotting eventually!
function plot_surfs(surfs)
    #expects an array of surface structs.

    display(surfs)

    p = plot(xlabel=L"r", ylabel=L"\phi", yguidefontrotation=0, left_margin=6Plots.mm, yguidefontsize=16, xguidefontsize=18, xtickfontsize=10, ytickfontsize=10, dpi=600, legendfontsize=10)

    for i in surfs

        s, t = convert_surf(i)

        plot!(t, s)

    end

    display(p)


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

    #display(length(surfs))
    ρsurfs = zeros(length(surfs))
    #so ρsurfs has to be in order
    #we can probably just have our inputs be in order, but 
    #we should also have some sorting here.
    #shouldn't be too hard.
    for (i, surf) in enumerate(surfs)
        #won't bother with the others yet.
        scosn[i, :, :] = surf.scos
        tsinn[i, :, :] = surf.tsin
        ρsurfs[i] = surf.ρ
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