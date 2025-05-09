"""
Struct storing the coeffiecients of the qfm surfaces.
Based on expansion 
r(ζ) = ∑ r_n^c cos(ζ) + r_n^s sin(ζ)
θ(ζ) = θ_0^c + bζ/a + ∑ θ_n^c cos(ζ) + θ_n^s sin(ζ)
"""
struct QFMSurfaceT
    rational :: Tuple{Int64, Int64} #(a, b)
    s :: Float64 #surface label
    rcos :: Array{Float64, 2}
    θsin :: Array{Float64, 2}
    rsin :: Array{Float64, 2}
    θcos :: Array{Float64, 2}
end

#name subject to change.
"""
Struct for storing the temporary values defining each surface used in coord_transform!
"""
struct TempSurfT
    rcos :: Array{Float64, 2} 
    θsin :: Array{Float64, 2}
    drcosds :: Array{Float64, 2}
    dθsinds :: Array{Float64, 2}
    d2rcosdsds :: Array{Float64, 2}
    d2θsindsds :: Array{Float64, 2}
    α :: Array{Float64, 2}
    cosα :: Array{Float64, 2}
    sinα :: Array{Float64, 2}
    function TempSurfT(M::Int64, N::Int64) 
        dim1 = M + 1
        dim2 = 2*N + 1
        new(zeros(dim1, dim2), zeros(dim1, dim2), zeros(dim1, dim2), zeros(dim1, dim2), zeros(dim1, dim2), zeros(dim1, dim2), zeros(dim1, dim2), zeros(dim1, dim2), zeros(dim1, dim2))
    end
end


"""
Struct storing the interpolations used to create surfaces in between the qfm surfaces found.
Derivatives are stored individually for extrapolation to work.
"""
struct SurfaceITPT
    M :: Int64
    N :: Int64
    rcos_itp :: Array{BSplineKit.SplineExtrapolations.SplineExtrapolation, 2}
    θsin_itp :: Array{BSplineKit.SplineExtrapolations.SplineExtrapolation, 2}
    drcos_itp :: Array{BSplineKit.SplineExtrapolations.SplineExtrapolation, 2}
    dθsin_itp :: Array{BSplineKit.SplineExtrapolations.SplineExtrapolation, 2}
    d2rcos_itp :: Array{BSplineKit.SplineExtrapolations.SplineExtrapolation, 2}
    d2θsin_itp :: Array{BSplineKit.SplineExtrapolations.SplineExtrapolation, 2}
end


"""
    construct_surfaces(rationals::Array{Tuple}, guesslist::Array{Float64}, prob::ProblemT; nfft=2::Int64, M=24::Int64, N=8::Int64)

Main function for computing the qfm surfaces. Iterates over the (a, b) rational inputs and creates a surface for each value.
"""
function construct_surfaces(rationals::Array{Tuple{Int64, Int64}}, guesslist::Array{Float64}, prob::ProblemT; nfft=2::Int64, M=24::Int64, N=8::Int64)


    met = MetT()
    B = BFieldT()
    surfaces = QFMSurfaceT[]

    for (i, rational) in enumerate(rationals)
        rcos, θsin, rsin, θcos = action(rational, prob, met, B, M, N, guesslist[i], nfft)

        s = rcos[1, 1] #surface label.
        push!(surfaces, QFMSurfaceT(rational, s, rcos, θsin, rsin, θcos))
        @printf("Found %d of %d surfaces.\n", i, length(rationals))
    end

    return surfaces
end



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
Creates a list of rationals, by getting all rationals between min and max, with a maximimum denominator of depth.
"""
function lowest_rationals(depth::Int64, min::Float64, max::Float64)

    rats = Tuple{Int64, Int64}[]
    for i in 1:depth
        j = 1
        while (j / i < max)
            if j/i > min
                a = gcd(j, i)
                if !((j, i) in rats)
                    if a==1
                        push!(rats, (j, i))
                    elseif !((Int(j/a), Int(i/a)) in rats)
                        push!(rats, (Int(j/a), Int(i/a)))
                    end
                end
            end
            j += 1
        end
    end
    return rats
end


"""

Creates a list of rationals using a farey tree.
"""
function farey_tree(N::Int64, q1::Int64, p1::Int64, q2::Int64, p2::Int64)

    #now N is a recursive depth.

    #think this is a stupid assertion.
    #won't help in weird cases,
    #also think it depends on the order.
    #@assert q2*p1 - q1*p2 == -1
    plist = [p1, p2]
    qlist = [q1, q2]

    farey_tree_recurs!(N, q1, p1, q2, p2, qlist, plist)

    #this enures the ratio is reduced.
    for i in eachindex(qlist)
        a = gcd(plist[i], qlist[i])
        if a != 1
            plist[i] ÷= a
            qlist[i] ÷= a
        end
    end

    #think we are going to assume that there are no double ups

    return qlist, plist

end

"""

Recursive part of farey_tree.
"""
function farey_tree_recurs!(N::Int64, q1::Int64, p1::Int64, q2::Int64, p2::Int64, qlist::Array{Int64}, plist::Array{Int64})

    if N != 0
        N = N-1
        pnew = p1 + p2
        qnew = q1 + q2
        push!(qlist, qnew)
        push!(plist, pnew)
        farey_tree_recurs!(N, qnew, pnew, q2, p2, qlist, plist)
        farey_tree_recurs!(N, q1, p1, qnew, pnew, qlist, plist)
    end


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
Guesses the best starting point for finding qfm surfaces based on the q_profile.
"""
function surface_guess(rationals::Array{Tuple}, q::Function)
    gl = Float64[]

    for i in rationals
        diff(r) = i[1]/i[2] - q(r)[1]
        sol = find_zero(diff, 0.5)
        push!(gl, sol)
    end
    return gl
end

"""
    compute_jac(prob::ProblemT, grids::FFFGridsT, surfs::Array{QFMSurfaceT})

Computes the Jacobain and Magnetic field in qfm coordinates to test the new values.
"""
function compute_jac(prob::ProblemT, grids::FFFGridsT, surfs::Array{QFMSurfaceT})

    #instantiate the grids into arrays. 
    #note that these grids are not actually the points used in construct, as we either use fourier exansion or gaussian weight points
    #but this should still give us a good idea, but care must be taken for r=0.
    rgrid, θgrid, ζgrid = inst_grids(grids)

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
    CT = CoordTsfmT()

    jac = zeros(grids.r.N, grids.θ.N, grids.ζ.N)
    djac = zeros(3, grids.r.N, grids.θ.N, grids.ζ.N)
    B = zeros(3, grids.r.N, grids.θ.N, grids.ζ.N)
    #for comparisons.
    #jac_tor = zeros(grids.r.N, grids.θ.N, grids.ζ.N)
    #djac_tor = zeros(3, grids.r.N, grids.θ.N, grids.ζ.N)
    #coords = zeros(3, grids.r.N, grids.θ.N, grids.ζ.N)

    #for (i, r) in enumerate(rvals), (j, θ) in enumerate(θvals), (k, ζ) in enumerate(ζvals)
    for (i, r) in enumerate(rgrid), (j, θ) in enumerate(θgrid), (k, ζ) in enumerate(ζgrid)
        coord_transform!(r, θ, ζ, CT, surf_itp, sd)
        toroidal_metric!(tor_met, CT.coords[1], CT.coords[2], CT.coords[3], prob.geo.R0)
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

##########################################################

#seems like just picking a rational surface near the boundary and running the og alg would be better?
#this doesn't seem to work in python either, perhaps we can just assume the perturbation is negligible at the boundaries and pick flat surfaces there?
#probably ok assuming the surfaces we do find are almost flat towards the edge
#NOTE: THIS DOES NOT WORK!
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

            prob.met(met, ρ, θarr[j, k], ζarr[1, j, k], prob.geo.R0)
            #dodge af. need to decide if this version of compute B will actually be used anywhere.
            compute_B!(B.B, met.J[1], prob.q, prob.isl, prob.isl2, ρ, θarr[j, k], ζarr[1, j, k])
            #zero chance this works.
            Bs[j, k] = B.B[1]
            Bθ[j, k] = B.B[2]
            Bζ[j, k] = B.B[3]
        end

        #this function is never going to work because of this step.
        #with our form of the magnetic field, this will always be a constant for a fixed ρ, then
        #the fourier transform will be nonsense. and the error goes to zero, because we get the same nonsense
        #answer twice in a row. Unsure how to fix this!.
        Bθ_o_Bζ = @. Bθ / Bζ

        cn, sn = rfft2D(Bθ_o_Bζ, M, N)

        #display(Bθ_o_Bζ)
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

        #display([erriota, errcn, errsn])
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
