
"""
    itp_mat!(surf_itp::SurfaceITPT, s::Float64)

Computes the matrix of interpolations.
"""
function itp_mat!(surf_itp::SurfaceITPT, sd::TempSurfT, s::Float64)

    dim1 = surf_itp.M + 1
    dim2 = 2*surf_itp.N + 1

    for i in 1:dim1, j in 1:dim2
        sd.rcos[i, j] = surf_itp.rcos_itp[i, j](s)
        sd.θsin[i, j] = surf_itp.θsin_itp[i, j](s)
        sd.drcosds[i, j] = surf_itp.drcos_itp[i, j](s)
        sd.dθsinds[i, j] = surf_itp.dθsin_itp[i, j](s)
        sd.d2rcosdsds[i, j] = surf_itp.d2rcos_itp[i, j](s)
        sd.d2θsindsds[i, j] = surf_itp.d2θsin_itp[i, j](s)
    end

end



"""
    create_surf_itp(surfs::Array{QFMSurfaceT})

Creates an a struct storing the interpolation function based on an array of qfm surfaces.
"""
function create_surf_itp(surfs::Array{QFMSurfaceT})

    dim1, dim2 = size(surfs[1].rcos)
    pqMpol = dim1 - 1
    pqNtor = (dim2 - 1) ÷ 2

    rcosn = zeros((length(surfs), dim1, dim2))
    θsinn = zeros((length(surfs), dim1, dim2))

    #sorts the surfaces by the ρ value, so that interpolation can work properly.
    s_surfs_unsort = zeros(length(surfs))

    for (i, surf) in enumerate(surfs)
        s_surfs_unsort[i] = surf.s
    end

    perm = sortperm(s_surfs_unsort)
    s_surfs = s_surfs_unsort[perm]
    surfs = surfs[perm]

    for (i, surf) in enumerate(surfs)
        rcosn[i, :, :] = surf.rcos
        θsinn[i, :, :] = surf.θsin
    end

    #arrays storing the interpolations.
    #we need individual arrays for the derivatives to ensure extrapolation works as intended.
    rcos_itp = Array{BSplineKit.SplineExtrapolations.SplineExtrapolation}(undef, dim1, dim2)
    drcos_itp = Array{BSplineKit.SplineExtrapolations.SplineExtrapolation}(undef, dim1, dim2)
    d2rcos_itp = Array{BSplineKit.SplineExtrapolations.SplineExtrapolation}(undef, dim1, dim2)
    θsin_itp = Array{BSplineKit.SplineExtrapolations.SplineExtrapolation}(undef, dim1, dim2)
    dθsin_itp = Array{BSplineKit.SplineExtrapolations.SplineExtrapolation}(undef, dim1, dim2)
    d2θsin_itp = Array{BSplineKit.SplineExtrapolations.SplineExtrapolation}(undef, dim1, dim2)


    for j in 1:dim2, i in 1:dim1

        rcos = interpolate(s_surfs, rcosn[:, i, j], BSplineOrder(5))
        θsin = interpolate(s_surfs, θsinn[:, i, j], BSplineOrder(5))
        rcos_itp[i, j] = extrapolate(rcos, Smooth())
        θsin_itp[i, j] = extrapolate(θsin, Smooth())
        #unclear exactly how well these will actually work outside the domain.
        #but will atleast stop linear algebra errors.
        drcos_itp[i, j] = extrapolate(Derivative(1)*rcos, Smooth())
        dθsin_itp[i, j] = extrapolate(Derivative(1)*θsin, Smooth())
        d2rcos_itp[i, j] = extrapolate(Derivative(2)*rcos, Smooth())
        d2θsin_itp[i, j] = extrapolate(Derivative(2)*θsin, Smooth())
    end

    #returns an empty struct of the correct size to store the intermediate interpolation data.
    return SurfaceITPT(pqMpol, pqNtor, rcos_itp, θsin_itp, drcos_itp, dθsin_itp, d2rcos_itp, d2θsin_itp), TempSurfT(pqMpol, pqNtor)

end


"""
Converts a surface into r, θ values.
"""
function convert_surf(surf::QFMSurfaceT)
    #takes the weird output form into a plotable form.


    Nϑ = 100
    ϑgrid = LinRange(0, 2*π, Nϑ)
    φ = 0.0

    #again, v stupid.
    rcos = surf.rcos
    rsin = surf.rsin
    θsin = surf.θsin
    θcos = surf.θcos

    α = zeros((Nϑ, size(rcos)[1], size(rcos)[2]))

    #this should be determinable from the array size.
    #if only we understood these arrays.
    #really good
    pqMpol, pqNtor = size(rcos)
    mlist = collect(range(0, pqMpol))

    collect(-pqNtor:0)
    nlist = [collect(0:pqNtor);collect(-pqNtor:-1)] 

    for i in 1:Nϑ
        for j in 1:size(rcos)[1]
            for k in 1:size(rcos)[2]
                α[i, j, k] = mlist[j] * ϑgrid[i] - nlist[k] * φ
            end
        end
    end

    cosα = cos.(α)
    sinα = sin.(α)

    r = zeros(Nϑ)
    θ = zeros(Nϑ)
    θ = collect(LinRange(0, 2π, Nϑ))

    for i in 1:Nϑ

        for j in 1:size(rcos)[1]
            for k in 1:size(rcos)[2]
                r[i] += rcos[j, k] * cosα[i, j, k]
                #θ[i] += ϑgrid[i] + θsin[j, k] * sinα[i, j, k]
                θ[i] += θsin[j, k] * sinα[i, j, k]
            end
        end
    end

    return r, θ

end


"""
    compute_jac(prob::ProblemT, grids::FFFGridsT, surfs::Array{QFMSurfaceT})

Computes the Jacobain and Magnetic field in qfm coordinates to test the new values.
"""
function compute_jac(prob::ProblemT, grids::FFFGridsT, surfs::Array{QFMSurfaceT})

    #instantiate the grids into arrays. 
    #note that these grids are not actually the points used in construct, as we either use fourier exansion or gaussian weight points
    #but this should still give us a good idea, but care must be taken for r=0.
    x1grid, x2grid, x3grid = inst_grids(grids)

    #initialise the two structs to store the metric and the magnetic field.
    tor_met = MetT()
    qfm_met = MetT()
    tor_B = BFieldT()
    qfm_B = BFieldT()

    #creates the interpolations for the surfaces.
    surf_itp, sd = create_surf_itp(surfs)

    #compute the gaussian qudrature points for finite elements.
    #ξr, wgr = MID.Construct.FastGaussQuadrature.gausslegendre(grids.r.gp) #same as python!
    #ξθ, wgθ = MID.Construct.FastGaussQuadrature.gausslegendre(grids.θ.gp)
    #ξζ, wgζ = MID.Construct.FastGaussQuadrature.gausslegendre(grids.ζ.gp)

    #struct for storing the intermediate data for the coordinate transform
    CT = CoordTransformT()

    jac = zeros(grids.x1.N, grids.x2.N, grids.x3.N)
    djac = zeros(3, grids.x1.N, grids.x2.N, grids.x3.N)
    B = zeros(3, grids.x1.N, grids.x2.N, grids.x3.N)
    #for comparisons.
    #jac_tor = zeros(grids.x1.N, grids.x2.N, grids.x3.N)
    #djac_tor = zeros(3, grids.x1.N, grids.x2.N, grids.x3.N)
    #coords = zeros(3, grids.x1.N, grids.x2.N, grids.x3.N)

    #for (i, r) in enumerate(rvals), (j, x2) in enumerate(x2vals), (k, x3) in enumerate(x3vals)
    for (i, x1) in enumerate(x1grid), (j, x2) in enumerate(x2grid), (k, x3) in enumerate(x3grid)
        coord_transform!(x1, x2, x3, CT, surf_itp, sd)
        prob.met(tor_met, CT.coords[1], CT.coords[2], CT.coords[3], prob.geo.R0)
        compute_B!(tor_B, tor_met, prob.q, prob.isls, CT.coords[1], CT.coords[2], CT.coords[3])
        met_transform!(tor_met, qfm_met, CT)
        B_transform!(tor_B, qfm_B, qfm_met, CT)

        jac[i, j, k] = qfm_met.J[1]
        djac[:, i, j, k] = qfm_met.dJ[:]
        #jac_tor[i, j, k] = tor_met.J[1]
        #djac_tor[:, i, j, k] = tor_met.dJ[:]
        #coords[:, i, j, k] = CT.coords[:]
        B[:, i, j, k] = qfm_B.B[:]

    end
    return B, jac, djac
end
