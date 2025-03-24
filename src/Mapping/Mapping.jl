

#handles coordinate mapping between island and non-island coords.
#think this at least mildy works,
#but this is awful
#requries a rewrite of island struct I think
#and redo of how island stuff is integrated into this code.
module Mapping


using Elliptic
using FFTW
#using Interpolations #one day we will just Hermite I think.
#need to swap interpolations to Bspline kit to match other case
using Printf
using JLD2


using MID.Structures
using MID.Geometry
using MID.Io
using MID.PostProcessing
using MID.Equilibrium


include("Interpolation.jl")
include("ToroidalToIsland.jl")
include("IslToToroidal.jl")

include("Continuum.jl")

export tor_to_isl_continuum
export isl_to_tor_continuum


export tor_to_isl
export isl_to_tor


#will need to make this more proper one day.
function isl_to_tor(tor_grids, ϕ, isl_grids, isl)

    #rgrid, θgrid, ζgrid = inst_grids(tor_grids)

    #κgrid, ᾱgrid, φgrid = inst_grids(isl_grids)
    rgrid = LinRange(0, tor_grids.κmax, tor_grids.Nκ)
    θgrid = LinRange(0, 2π, tor_grids.Nᾱ)
    ζgrid = LinRange(0, 2π, tor_grids.Nφ)

    #so we actually have periodic grids.
    #κgrid = LinRange(0, isl_grids.κmax, isl_grids.Nκ)
    #this is extremely confusing because the grids are designed
    #for the other way...
    #so all vars are kinda switched.
    κgrid = LinRange(0, 1.0, isl_grids.r.N)
    ᾱgrid = LinRange(0, 2π, isl_grids.θ.N)
    φgrid = LinRange(0, 2π, isl_grids.ζ.N)

    #will ignore periodicty for now!

    isl_itp = interpolate((κgrid, ᾱgrid, φgrid), ϕ, (Gridded(Linear()), Gridded(Linear(Periodic())), Gridded(Linear(Periodic()))));

    isl_ext = extrapolate(isl_itp, Periodic());


    ϕ_tor = Array{ComplexF64}(undef, length(rgrid), length(θgrid), length(ζgrid));


    for (i, r) in enumerate(rgrid), (j, θ) in enumerate(θgrid), (k, ζ) in enumerate(ζgrid)

        κ, ᾱ, φ = coords_tor_to_isl(r, θ, ζ, isl)

        if κ==0 && ᾱ==0 && φ==0
            #for some reason this case is being mapped to 2?
            ϕ_tor[i, j, k] = 0
        else

            ϕ_tor[i, j, k] = isl_ext(κ, ᾱ, φ)
        end
    end

    return ϕ_tor, fft(ϕ_tor, [2, 3])



end

#lets muck around with the inside vs outside domain etc.
#best method may be to restrict ϕ to only focus on a single island??? 
#that will be quite difficult though.
function new_tor_to_isl(mapgrids::MapGridsT, ϕ, grids::FFFGridsT, isl::IslandT, sign=1)

    #this uses inbuilt interpolation, which I think helps with domain issues
    #also we use Axel's formulation.

    #note that the ϕ to be mapped must be fff and contain the deriv parts.
    #may implement others some day.

    #first we define the equivalent island grids

    rgrid, θgrid, ζgrid = inst_grids(grids)

    #I think this is bad, but hopefully we won't need to do this once we get Hermite working.
    ϕ_p = periodify(ϕ, grids.r.N, grids.θ.N, grids.ζ.N)

    θgrid_p = LinRange(0, 2π, grids.θ.N + 1)
    ζgrid_p = LinRange(0, 2π, grids.ζ.N + 1)

    #ideally swap to hermite interpolation at some stage.
    tor_itp = interpolate((rgrid, θgrid_p, ζgrid_p), ϕ_p, (Gridded(Linear()), Gridded(Linear(Periodic())), Gridded(Linear(Periodic()))));

    tor_ext = extrapolate(tor_itp, Periodic());

    #not sure how we are going to deal with κ max. may just need to be a parameter.
    #2 is completly arbitrary lol.
    κgrid, ᾱgrid, φgrid = inst_grids(mapgrids)

    #ᾱgrid = LinRange(0, π, mapgrids.Nᾱ) #test

    κgrid_in = sqrt.(κgrid[κgrid .< 1])
    κgrid_out = sqrt.(κgrid[κgrid .>= 1])

    ϕ_islp = Array{ComplexF64}(undef, mapgrids.Nκ, mapgrids.Nᾱ, mapgrids.Nφ);

    ϕ_islm = Array{ComplexF64}(undef, mapgrids.Nκ, mapgrids.Nᾱ, mapgrids.Nφ);

    ϕ_isl_in = Array{ComplexF64}(undef, length(κgrid_in), mapgrids.Nᾱ, mapgrids.Nφ);
    ϕ_isl_outp = Array{ComplexF64}(undef, length(κgrid_out), mapgrids.Nᾱ, mapgrids.Nφ);
    ϕ_isl_outm = Array{ComplexF64}(undef, length(κgrid_out), mapgrids.Nᾱ, mapgrids.Nφ);


    #do we need to map inside and outside separatly??? Sounds like a terrible solution.

    #we probbaly need to redefine α to match the definition of our islands.

    #iterate through each point on our island grid
    for (i, κ) in enumerate(κgrid_in), (j, ᾱ) in enumerate(ᾱgrid), (k, φ) in enumerate(φgrid)

        #find the equivalent coordinate in toroidal coordinates.
        #fkn worst name so far.
        r, θ, ζ = coords_isl_to_tor(κ, ᾱ, φ, isl, sign)

        #now we interpolate this coordinate using our functions.

        ϕ_isl_in[i, j, k] = tor_ext(r, θ, ζ)
        ϕ_islp[i, j, k] = tor_ext(r, θ, ζ)
        ϕ_islm[i, j, k] = tor_ext(r, θ, ζ)
        

    end

    ᾱgrid = LinRange(0, π, mapgrids.Nᾱ) #test

    for (i, κ) in enumerate(κgrid_out), (j, ᾱ) in enumerate(ᾱgrid), (k, φ) in enumerate(φgrid)

        #find the equivalent coordinate in toroidal coordinates.
        #fkn worst name so far.
        r, θ, ζ = coords_isl_to_tor(κ, ᾱ, φ, isl, +1)

        #now we interpolate this coordinate using our functions.

        ϕ_isl_outp[i, j, k] = tor_ext(r, θ, ζ)
        ϕ_islp[i+length(κgrid_in), j, k] = tor_ext(r, θ, ζ)
        

        r, θ, ζ = coords_isl_to_tor(κ, ᾱ, φ, isl, -1)

        #now we interpolate this coordinate using our functions.
        ϕ_isl_outm[i, j, k] = tor_ext(r, θ, ζ)
        ϕ_islm[i+length(κgrid_in), j, k] = tor_ext(r, θ, ζ)
        
        

    end
    #fft(ϕ_isl, [2, 3])
    return ϕ_isl_in, ϕ_isl_outp, ϕ_isl_outm, ϕ_islp, ϕ_islm
    

end


#unmodified version...
function tor_to_isl(mapgrids::MapGridsT, ϕ, grids::FFFGridsT, isl::IslandT, sign=1)

    #this uses inbuilt interpolation, which I think helps with domain issues
    #also we use Axel's formulation.

    #note that the ϕ to be mapped must be fff and contain the deriv parts.
    #may implement others some day.

    #first we define the equivalent island grids

    rgrid, θgrid, ζgrid = inst_grids(grids)

    #I think this is bad, but hopefully we won't need to do this once we get Hermite working.
    ϕ_p = periodify(ϕ, grids.r.N, grids.θ.N, grids.ζ.N)

    θgrid_p = LinRange(0, 2π, grids.θ.N + 1)
    ζgrid_p = LinRange(0, 2π, grids.ζ.N + 1)

    #ideally swap to hermite interpolation at some stage.
    tor_itp = interpolate((rgrid, θgrid_p, ζgrid_p), ϕ_p, (Gridded(Linear()), Gridded(Linear(Periodic())), Gridded(Linear(Periodic()))));

    tor_ext = extrapolate(tor_itp, Periodic());

    #not sure how we are going to deal with κ max. may just need to be a parameter.
    #2 is completly arbitrary lol.
    κgrid, ᾱgrid, φgrid = inst_grids(mapgrids)

    ϕ_isl = Array{ComplexF64}(undef, mapgrids.Nκ, mapgrids.Nᾱ, mapgrids.Nφ);


    #do we need to map inside and outside separatly??? Sounds like a terrible solution.

    #we probbaly need to redefine α to match the definition of our islands.

    #iterate through each point on our island grid
    for (i, κ) in enumerate(κgrid), (j, ᾱ) in enumerate(ᾱgrid), (k, φ) in enumerate(φgrid)

        #find the equivalent coordinate in toroidal coordinates.
        #fkn worst name so far.
        r, θ, ζ = coords_isl_to_tor(κ, ᾱ, φ, isl, sign)

        #now we interpolate this coordinate using our functions.

        ϕ_isl[i, j, k] = tor_ext(r, θ, ζ)
        
        

    end
    
    return ϕ_isl, fft(ϕ_isl, [2, 3])
    

end


#should use grids.
#will be annoying af if we end up doing the Hermite version of this...
#unclear exactly how important this really is.
function periodify(ϕ, Nr, Nθ, Nζ)

    ϕ_p = zeros(ComplexF64, Nr, Nθ+1, Nζ+1)

    ϕ_p[:, 1:end-1, 1:end-1] = ϕ
    ϕ_p[:, end, end] = ϕ[:, 1, 1]

    return ϕ_p

end



function tor_to_isl_Axel(Nκ, Nβs, Nφ, ϕ, grids::FFFGridsT, isl::IslandT)

    #this uses inbuilt interpolation, which I think helps with domain issues
    #also we use Axel's formulation.

    #note that the ϕ to be mapped must be fff and contain the deriv parts.
    #may implement others some day.

    #first we define the equivalent island grids

    rgrid, θgrid, ζgrid = instantiate_grids(grids)


    tor_itp = interpolate((rgrid, θgrid, ζgrid), ϕ, (Gridded(Linear()), Gridded(Linear(Periodic())), Gridded(Linear(Periodic()))));

    tor_ext = extrapolate(tor_itp, Periodic());

    #not sure how we are going to deal with κ max. may just need to be a parameter.
    #2 is completly arbitrary lol.
    κgrid = LinRange(0, 2, Nκ);
    #κgrid2 = LinRange(0, 1.5, Nκ2);
    βsgrid = LinRange(-π, π, Nβs+1)[1:end-1];

    #changing this to -pi to pi didn't really change anything, just scaled some stuff.
    #which means this domain does have an effect on ft... not ideal
    φgrid = LinRange(0, 2π, Nφ+1)[1:end-1];

    ϕ_isl = zeros(ComplexF64, Nκ, Nβs, Nφ);

    #iterate through each point on our island grid
    for (i, κ) in enumerate(κgrid), (j, βs) in enumerate(βsgrid), (k, φ) in enumerate(φgrid)

        #find the equivalent coordinate in toroidal coordinates.
        #fkn worst name so far.
        r, θ, ζ = coords_isl_to_tor(κ, βs, φ, isl)

        #now we interpolate this coordinate using our functions.

        ϕ_isl[i, j, k] = tor_ext(r, θ, ζ)
        
        

    end

    #perhaps we should fft this???
    
    return ϕ_isl, fft(ϕ_isl, [2, 3])
    

end


#this maps a toroidal function into an island function.
#going to avoid the hermite stuff now to try and reduce error areas.
function tor_to_isl_hermite(Nκ, Nᾱ, Nφ, ϕ, grids::FFFGridsT, isl::IslandT)

    #note that the ϕ to be mapped must be fff and contain the deriv parts.
    #may implement others some day.

    #first we define the equivalent island grids

    #not sure how we are going to deal with κ max. may just need to be a parameter.
    #2 is completly arbitrary lol.
    κgrid = LinRange(0, 2, Nκ);
    #κgrid2 = LinRange(0, 1.5, Nκ2);
    ᾱgrid = LinRange(0, 2π, Nᾱ+1)[1:end-1];
    φgrid = LinRange(0, 2π, Nφ+1)[1:end-1];

    ϕ_isl = zeros(ComplexF64, Nκ, Nᾱ, Nφ);

    #iterate through each point on our island grid
    for (i, κ) in enumerate(κgrid), (j, ᾱ) in enumerate(ᾱgrid), (k, φ) in enumerate(φgrid)

        #find the equivalent coordinate in toroidal coordinates.
        #fkn worst name so far.
        r, θ, ζ = coords_isl_to_tor(κ, ᾱ, φ, isl)

        #now we interpolate this coordinate using our functions.

        ϕ_isl[i, j, k] = hermite_interpolation(r, mod(θ, 2π), mod(ζ, 2π), ϕ, grids)

    end

    #perhaps we should fft this???
    
    return ϕ_isl, fft(ϕ_isl, [2, 3])
    

end

end
