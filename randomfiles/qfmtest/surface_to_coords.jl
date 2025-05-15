
#once surfaces have actually been created, this file will make them actually useful
using Interpolations

#new interpolations package to test
#this should be the same as scipy.
using Dierckx


#struct for storing the surface interpolations.
#class based structure sure is nice for this application.
#this will almost certainly require more.
struct surface_itp
    scos_itp
    tsin_itp
end



#these are the 4 things the define a surface I guess??
#we will want to rewrite these functions to be a lot clearer.
#we probably want an mlist and nlist creating method/storage.
#unsure if we in general want to consider only a single coord input
#or if we want to do everthing at once???
#perhaps we assume single coord at a time.
#may be worth implementing a coord struct at some stage, may make our code cleaner, then we can probably use r, θ, ζ or κ, α etc without much difference
#probably useful here as we will be introducing new coords again...
#this is a stupid af function, but perhaps it is worth doing this way so we can see what each function will need for when we restructrue.

#this needs to do more than initially thought.
#so other than interpolations not working
#and being awful
#and probably has typos
#this function does what is needed
#i.e. it return the correct shape at least!
function transform_coords(r, θ, ζ, surf_itp)

    #this is rapidly becoming a complete disaster of a function.
    #we probably want to split different parts up if possible
    #also this is in the running for worst function ever written.

    #obvs not ideal...
    pqMpol = 24
    pqNtor = 8
    mlist = collect(range(0, pqMpol))

    collect(-pqNtor:0)
    nlist = [collect(0:pqNtor);collect(-pqNtor:-1)] 

    #this is going to be v similar to convert_surf,
    #just without all the points matching identically.

    #start by assuming only a single input coord.
    α = zeros((length(mlist), length(nlist)))

    for i in 1:1:length(mlist)
        for j in 1:1:length(nlist)
            α[i, j] = mlist[i] * θ - nlist[j] * ζ
        end
    end

    cosα = cos.(α)
    sinα = sin.(α)

    #currently restricted for r values between our outermost flux surfaces.
    scos = surf_itp.scos_itp(r)
    tsin = surf_itp.tsin_itp(r)

    #this may not work
    dscosds = only(Interpolations.gradient(surf_itp.scos_itp, r))
    dtsinds = only(Interpolations.gradient(surf_itp.tsin_itp, r))

    #this almost certainly won't work.
    #probably should try the Dierckx.jl package, which wraps the same thing scipy interpolations wraps.
    #should also be able to do cubic stuff if needed.
    #wot the fek even is s. So unclear.
    #calling this is giving a stack overflow error lol.
    #unsure what is going on here lol.
    #may mean we need to use the other package.
    #thiink because this is using linear shit anyway,
    #this won't work.
    #d2scosds2 = only(Interpolations.hessian(surf_itp.scos_itp, r))
    #d2tsinds2 = only(Interpolations.hessian(surf_itp.tsin_itp, r))

    d2scosds2 = zeros(25, 17)
    d2tsinds2 = zeros(25, 17)

    msinα = zeros(size(sinα))
    mcosα = zeros(size(cosα))
    nsinα = zeros(size(sinα))
    ncosα = zeros(size(cosα))

    #second derivatives.
    m2sinα = zeros(size(sinα))
    m2cosα = zeros(size(cosα))
    n2sinα = zeros(size(sinα))
    n2cosα = zeros(size(cosα))
    mnsinα = zeros(size(sinα))
    mncosα = zeros(size(cosα))
    
    

    #we may need to use a different package to get the second derivative to work
    s = 0
    t = θ

    for i in 1:1:length(mlist)
        for j in 1:1:length(nlist)
            s += scos[i, j] * cosα[i, j]
            t += tsin[i, j] * sinα[i, j]
            #these are a complete waste of memory,
            #we are iterating anyway.
            #no point doing this.
            msinα[i, j] = mlist[i] * sinα[i, j]
            mcosα[i, j] = mlist[i] * cosα[i, j]
            nsinα[i, j] = nlist[j] * sinα[i, j]
            ncosα[i, j] = nlist[j] * cosα[i, j]

            m2sinα[i, j] = mlist[i]^2 * sinα[i, j]
            m2cosα[i, j] = mlist[i]^2 * cosα[i, j]
            n2sinα[i, j] = nlist[j]^2 * sinα[i, j]
            n2cosα[i, j] = nlist[j]^2 * cosα[i, j]
            mnsinα[i, j] = mlist[i] * nlist[j] * sinα[i, j]
            mncosα[i, j] = mlist[i] * nlist[j] * cosα[i, j]


        end
    end

    dsdr = 0.0
    dtdr = 0.0
    dsdt = 0.0
    #note this starts at 1!!!
    dtdt = 1.0 #srs wot the fek.

    dsdz = 0.0
    dtdz = 0.0

    d2sdrdr = 0.0
    d2sdrdt = 0.0
    d2sdrdz = 0.0
    d2sdtdt = 0.0
    d2sdtdz = 0.0
    d2sdzdz = 0.0
    
    d2tdrdr = 0.0
    d2tdrdt = 0.0
    d2tdrdz = 0.0
    d2tdtdt = 0.0
    d2tdtdz = 0.0
    d2tdzdz = 0.0


    for i in 1:1:length(mlist)
        for j in 1:1:length(nlist)
            #this is a disgrace.
            dsdr += dscosds[i, j] * cosα[i, j]
            dtdr += dtsinds[i, j] * sinα[i, j]
            #given we are already looping through this, we do not need
            #to actually store msin etc. #TODO.
            dsdt -= scos[i, j] * msinα[i, j]
            dtdt += tsin[i, j] * mcosα[i, j]
            dsdz -= ssin[i, j] * ncosα[i, j]
            dtdz += tcos[i, j] * nsinα[i, j]


            #now we just have to do the second derivs.....
            #we will need to check all of these. This is a disaster.
            d2sdrdr += d2scosds2[i, j] * cosα[i, j]

            d2sdrdt -= dscosds[i, j] * msinα[i, j]

            #plus here means we are using a different convention for m+n in fourier here. not ideal.
            d2sdrdz += dscosds[i, j] * nsinα[i, j]

            d2sdtdt -= scos[i, j] * m2cosα[i, j]

            d2sdtdz += scos[i, j] * mncosα[i, j]

            d2sdzdz -= scos[i, j] * n2cosα[i, j]

            #now theta

            d2tdrdr += d2tsinds2[i, j] * sinα[i, j]

            d2tdrdt -= dtsinds[i, j] * mcosα[i, j]

            
            d2tdrdz -= dtsinds[i, j] * ncosα[i, j]

            d2tdtdt -= tsin[i, j] * m2sinα[i, j]

            d2tdtdz += tsin[i, j] * mnsinα[i, j]

            d2tdzdz -= tsin[i, j] * n2sinα[i, j]


        end
    end
    #other derivatives are computed via fourier derives.
    

    #so to be clear, t is transforming to t.
    #but r is troansforming to s. fkn stupid
    #once we have a better grasp of this we need to rename everything.
    J = dsdr * dtdt - dtdr * dsdt

    #properly have no idea what equation Zhisong is doing for this.
    dJ = zeros(3)

    #maybe correct, who knows.
    dJ[1] = dsdr * d2tdrdt + d2sdrdr * dtdt - (dtdr * d2sdrdt + d2tdrdr * dsdt)
    dJ[2] = dsdr * d2tdtdt + d2sdrdt * dtdt - (dtdr * d2sdtdt + d2tdrdt * dsdt)
    dJ[3] = dsdr * d2tdtdz + d2sdrdz * dtdt - (dtdr * d2sdtdz + d2tdrdz * dsdt)

    #order of this is probably cooked.
    JM = [dsdr dsdt dsdz ; dtdr dtdt dtdz ; 0.0 0.0 1.0]

    # this is cooked.
    #needs to be a 3x3x3 matrix, unsure how best to enter this tbh.
    #dJM = 


    return s, t, ζ, J, dJ, JM
end




#this is almost certainly not good enough!
function create_interpolation(surfs)
    #surfs is an array of the surface struct, which needs some serious fixing.
    #so we need some rhosurfs for this???
    #hard copied from python, this are found during qfm construction I guess.
    #note that python also has the 'edges' problemo for another day I think.

    #hint that our data storage method is not ideal
    scosn = zeros((length(surfs), size(surfs[1].scos)[1], size(surfs[1].scos)[2]))
    tsinn = zeros((length(surfs), size(surfs[1].scos)[1], size(surfs[1].scos)[2]))

    display(length(surfs))
    for (i, surf) in enumerate(surfs)
        #won't bother with the others yet.
        scosn[i, :, :] = surf.scos
        tsinn[i, :, :] = surf.tsin
    end

    ρsurfs = [0.6, 0.6111111111, 0.61538462, 0.61904762, 0.625]

    pqMpol = 24
    pqNtor = 8
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
    scosn = [surfs[i].scos for i in 1:length(ρsurfs)]
    tsinn = [surfs[i].tsin for i in 1:length(ρsurfs)]

    #this is satisfactory for now, but this function needs work
    #and we will need to compare more direclty with Zhisong
    #to make sure this form of interpolation is good enough!
    scos_itp = interpolate((ρsurfs, ), scosn, Gridded(Linear()))
    tsin_itp = interpolate((ρsurfs, ), tsinn, Gridded(Linear()))


    #display(scos_itp(0.62))
    #ok, so to get this to work.
    #we want to interpolate θ, ζ (or their fourier values or whatever)
    #based on a single ρ value.
    #so we need to rewrite scosn to be an array of length ρ length
    #and each point contains the θ and ζ values.
    #we may need to flatten the scosn array in θ and ζ
    #scos_itp = interpolate((Base.OneTo(5), ρsurfs, ρsurfs), scosn, (NoInterp(), Gridded(Linear()), Gridded(Linear())))
    #, BSpline())
    #tsin_itp = interpolate(tsinn, BSpline())

    return surface_itp(scos_itp, tsin_itp)




end

#we will need to test this with know coord transformation I think, so that we can be sure it is working
#perhaps compute some value in spherical and cylindrical
#should just require that we modify surf_itp. hopefully not to hard.
function qfm_geometry(r, θ, ζ, met, qfm_met, surf_itp)
    #this will also assume a single coord at a time.
    #input is the og coordinates
    #which we are mapping to the appropriate flux coordinate.
    #ideally the flux coords will have differnet labels
    #unsure what they should be.

    #we will denote the flux coords as.
    #(ρ, ϑ, φ) 

    R0 = 10.0 #TODO

    #obvs will want this to be a choice.
    toroidal_met!(met, r, θ, ζ, R0)

    #coord J is like the jacobian of the transformation.
    #this is going to be real confusing real quick
    #coords_JM is the full jacobian matrix of the transformation
    #i.e. the 3x3 matrix with all the derivs, needed to compute the metric.
    #need some better names. Having multiple jacobians that are actually different is no good.
    #I am not sure this should call transform_coords
    #as we don't actually use the new coords here, just the extra stuff
    #may be better to call transform coords elsewhere and pass the extra stuff in.
    #we will need to work on the global structure a lot.
    ρ, ϑ, φ, coord_J, coord_dJ, coords_JM = transform_coords(r, θ, ζ, surf_itp)

    qfm_met.J = met.J * coord_J
    qfm_met.dJ = met.dJ .* coord_J + met.J .* coord_dJ

    #I am not sure if we will be able to compute dgl and dgu...
    #that may require the inverse? of coords_JM. We might actually have to think about this coordinate transformation damn.
    #need to check these formulas etc.
    #this should be relativly straightforward to do an example from cylindrical to spherical.
    qfm_met.gl = coords_JM' * met.gl * coords_JM

    qfm_met.gu = inverse(qfm_met.gl)

    #we will need to work out how to do the derivatives ourselves RIP.
    #dgl shouldn't be to bad
    #but dgu may be v difficult.
    #maybe the best way to do it, is to write the analytical terms of gu in terms of gl, then take the derivative?
    #will be awful but doable.
    #unless there is a better relationship somewhere.

    return qfm_met
    
end