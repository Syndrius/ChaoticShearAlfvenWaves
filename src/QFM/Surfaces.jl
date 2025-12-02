
"""
    itp_mat!(surf_itp::SurfaceITPT, sd::TempSurftT, s::Float64)

Computes the matrix of interpolations.
"""
function itp_mat!(surf_itp::SurfaceITPT, sd::TempSurfT, s::Float64)

    #dimension is based on Fourier expansion, see CSAWCantori.
    dim1 = surf_itp.M + 1
    dim2 = 2*surf_itp.N + 1

    for j in 1:dim2, i in 1:dim1
        sd.ψcos[i, j] = surf_itp.ψcos_itp[i, j](s)
        sd.θsin[i, j] = surf_itp.θsin_itp[i, j](s)
        sd.dψcosds[i, j] = surf_itp.dψcos_itp[i, j](s)
        sd.dθsinds[i, j] = surf_itp.dθsin_itp[i, j](s)
        sd.d2ψcosdsds[i, j] = surf_itp.d2ψcos_itp[i, j](s)
        sd.d2θsindsds[i, j] = surf_itp.d2θsin_itp[i, j](s)
    end

end



"""
    create_surf_itp(surfs::Array{QFMSurfaceT})

Creates an a struct storing the interpolation function based on an array of qfm surfaces.
"""
function create_surf_itp(surfs::Array{QFMSurfaceT})

    dim1, dim2 = size(surfs[1].ψcos)
    pqMpol = dim1 - 1
    pqNtor = (dim2 - 1) ÷ 2

    ψcosn = zeros((length(surfs), dim1, dim2))
    θsinn = zeros((length(surfs), dim1, dim2))

    #sorts the surfaces by the s value, so that interpolation can work properly.
    s_surfs_unsort = zeros(length(surfs))

    for (i, surf) in enumerate(surfs)
        s_surfs_unsort[i] = surf.s
    end

    perm = sortperm(s_surfs_unsort)
    s_surfs = s_surfs_unsort[perm]
    surfs = surfs[perm]

    for (i, surf) in enumerate(surfs)
        ψcosn[i, :, :] = surf.ψcos
        θsinn[i, :, :] = surf.θsin
    end

    #arrays storing the interpolations.
    #we need individual arrays for the derivatives to ensure extrapolation works as intended.
    ψcos_itp = Array{BSplineKit.SplineExtrapolations.SplineExtrapolation}(undef, dim1, dim2)
    dψcos_itp = Array{BSplineKit.SplineExtrapolations.SplineExtrapolation}(undef, dim1, dim2)
    d2ψcos_itp = Array{BSplineKit.SplineExtrapolations.SplineExtrapolation}(undef, dim1, dim2)
    θsin_itp = Array{BSplineKit.SplineExtrapolations.SplineExtrapolation}(undef, dim1, dim2)
    dθsin_itp = Array{BSplineKit.SplineExtrapolations.SplineExtrapolation}(undef, dim1, dim2)
    d2θsin_itp = Array{BSplineKit.SplineExtrapolations.SplineExtrapolation}(undef, dim1, dim2)


    for j in 1:dim2, i in 1:dim1

        ψcos = interpolate(s_surfs, ψcosn[:, i, j], BSplineOrder(5))
        θsin = interpolate(s_surfs, θsinn[:, i, j], BSplineOrder(5))
        ψcos_itp[i, j] = extrapolate(ψcos, Smooth())
        θsin_itp[i, j] = extrapolate(θsin, Smooth())
        #unclear exactly how well these will actually work outside the domain.
        #but will atleast stop linear algebra errors.
        dψcos_itp[i, j] = extrapolate(Derivative(1)*ψcos, Smooth())
        dθsin_itp[i, j] = extrapolate(Derivative(1)*θsin, Smooth())
        d2ψcos_itp[i, j] = extrapolate(Derivative(2)*ψcos, Smooth())
        d2θsin_itp[i, j] = extrapolate(Derivative(2)*θsin, Smooth())
    end

    #returns an empty struct of the correct size to store the intermediate interpolation data.
    return SurfaceITPT(pqMpol, pqNtor, ψcos_itp, θsin_itp, dψcos_itp, dθsin_itp, d2ψcos_itp, d2θsin_itp), TempSurfT(pqMpol, pqNtor)

end


"""
    convert_surf(surf::QFMSurfaceT)

Converts a surface into ψ, θ values for plotting.
"""
function convert_surf(surf::QFMSurfaceT)

    Nϑ = 100
    ϑgrid = LinRange(0, 2*π, Nϑ)
    ζ = 0.0

    ψcos = surf.ψcos
    ψsin = surf.ψsin
    θsin = surf.θsin
    θcos = surf.θcos

    α = zeros((Nϑ, size(ψcos)[1], size(ψcos)[2]))

    pqMpol, pqNtor = size(ψcos)
    mlist = collect(range(0, pqMpol))

    collect(-pqNtor:0)
    nlist = [collect(0:pqNtor);collect(-pqNtor:-1)] 

    for i in 1:Nϑ
        for j in 1:size(ψcos)[1]
            for k in 1:size(ψcos)[2]
                α[i, j, k] = mlist[j] * ϑgrid[i] - nlist[k] * ζ
            end
        end
    end

    cosα = cos.(α)
    sinα = sin.(α)

    ψ = zeros(Nϑ)
    θ = collect(LinRange(0, 2π, Nϑ))

    for i in 1:Nϑ

        for j in 1:size(ψcos)[1]
            for k in 1:size(ψcos)[2]
                ψ[i] += ψcos[j, k] * cosα[i, j, k]
                θ[i] += θsin[j, k] * sinα[i, j, k]
            end
        end
    end

    return ψ, θ

end


"""
    compute_jac(prob::ProblemT, grids::FFFGridsT, surfs::Array{QFMSurfaceT})

Computes the Jacobian and Magnetic field in qfm coordinates to test the new values.
This is used for deciding on QFM surfaces.
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
        prob.geo.met(tor_met, CT.coords[1], CT.coords[2], CT.coords[3], prob.geo.R0)
        compute_B!(tor_B, tor_met, prob.fields.q, prob.fields.isls, CT.coords[1], CT.coords[2], CT.coords[3])
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
