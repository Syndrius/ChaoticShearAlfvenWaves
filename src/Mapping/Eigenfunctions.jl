

#this is the same for each case.
function efunc_map!(ϕ_map::Array{ComplexF64, 3}, N1::Int64, N2::Int64, N3::Int64, ϕ::Array{ComplexF64, 4}, x1g::AbstractArray{Float64}, x2g::AbstractArray{Float64}, x3g::AbstractArray{Float64}, coord_map::Array{Tuple{Float64, Float64, Float64}})

    for i in 1:N1, j in 1:N2, k in 1:N3
        ϕ_map[i, j, k] = hermite_interpolation(coord_map[i, j, k]..., ϕ, x1g, x2g, x3g)
    end
end

#should put the coord map construction in here!


############################### #all old now!
function qfm_efunc_to_isl!(ϕ_isl::Array{ComplexF64, 3}, κgrid::AbstractArray{Float64}, ᾱgrid::AbstractArray{Float64}, τgrid::AbstractArray{Float64}, ϕ_qfm::Array{ComplexF64, 3}, sgrid::AbstractArray{Float64}, ϑgrid::AbstractArray{Float64}, φgrid::AbstractArray{Float64}, isl::IslandT, CT::CoordTransformT, surf_itp::SurfaceITPT, sd::TempSurfT) 


    qfm_itp = interpolate((sgrid, ϑgrid, φgrid), ϕ_qfm, (Gridded(Linear()), Gridded(Linear(Periodic())), Gridded(Linear(Periodic()))));

    qfm_ext = extrapolate(qfm_itp, Periodic());

    #This will be slow af unfort
    #imagine if we still used Hermite interpolation!
    for (i, κ) in enumerate(κgrid), (j, ᾱ) in enumerate(ᾱgrid), (k, τ) in enumerate(τgrid)

        #first change to tor coords
        r, θ, ζ = isl_in_coords_to_tor(κ, ᾱ, τ, isl)
        #then change to qfm
        s, ϑ, φ = tor_coords_to_qfm(r, θ, ζ, CT, surf_itp, sd)

        ϕ_isl[i, j, k] = qfm_ext(s, ϑ, φ)

    end


end


function qfm_efunc_to_tor!(ϕ_tor::Array{ComplexF64, 3}, Nr::Int64, Nθ::Int64, Nζ::Int64, ϕ_qfm::Array{ComplexF64, 4}, sgrid::AbstractArray{Float64}, ϑgrid::AbstractArray{Float64}, φgrid::AbstractArray{Float64}, coord_map::Array{Tuple{Float64, Float64, Float64}})


    #qfm_itp = interpolate((sgrid, ϑgrid, φgrid), ϕ_qfm, (Gridded(Linear()), Gridded(Linear(Periodic())), Gridded(Linear(Periodic()))));

    #qfm_ext = extrapolate(qfm_itp, Periodic());

    #we may be able to optimise this by noting that ζ=φ, so we can kind of skip that loop
    #same for the island case.
    #although, in both cases, the new poloidal angle is a function of the toroidal angle.
    #can probably change this loop to be without the rgrid etc now.
    #would also help the function args be much nicer
    for i in 1:Nr, j in 1:Nθ, k in 1:Nζ
        #this is unfort a bit fked.
        #assuming this work ok though, this should be fine.
        #s, ϑ, φ = tor_coords_to_qfm(r, θ, ζ, CT, surf_itp, sd)

        #ϕ_tor[i, j, k] = qfm_ext(s, ϑ, φ)
        #ϕ_tor[i, j, k] = qfm_ext(coord_map[i, j, k]...)
        ϕ_tor[i, j, k] = hermite_interpolation(coord_map[i, j, k]..., ϕ_qfm, sgrid, ϑgrid, φgrid)
        #ϕ_isl[i, j, k] = hermite_interpolation(r, θ, ζ, ϕ_tor, rgrid, θgrid, ζgrid)

    end


end

function qfm_efunc_to_tor!(ϕ_tor::Array{ComplexF64, 3}, rgrid::AbstractArray{Float64}, θgrid::AbstractArray{Float64}, ζgrid::AbstractArray{Float64}, ϕ_qfm::Array{ComplexF64, 3}, sgrid::AbstractArray{Float64}, ϑgrid::AbstractArray{Float64}, φgrid::AbstractArray{Float64}, coord_map::Array{Vector{Float64}})


    qfm_itp = interpolate((sgrid, ϑgrid, φgrid), ϕ_qfm, (Gridded(Linear()), Gridded(Linear(Periodic())), Gridded(Linear(Periodic()))));

    qfm_ext = extrapolate(qfm_itp, Periodic());

    #we may be able to optimise this by noting that ζ=φ, so we can kind of skip that loop
    #same for the island case.
    #although, in both cases, the new poloidal angle is a function of the toroidal angle.
    for (i, r) in enumerate(rgrid), (j, θ) in enumerate(θgrid), (k, ζ) in enumerate(ζgrid)
        #this is unfort a bit fked.
        #assuming this work ok though, this should be fine.
        #s, ϑ, φ = tor_coords_to_qfm(r, θ, ζ, CT, surf_itp, sd)

        #ϕ_tor[i, j, k] = qfm_ext(s, ϑ, φ)
        ϕ_tor[i, j, k] = qfm_ext(coord_map[i, j, k]...)

    end


end

#grids probbaly need to be fff for now, mainly because the fft extends the grids.
#island also needs to be instantiated!
#this function takes an eigenfunction computed using island coordinates, and maps it to toroidal coordinates.
function isl_efunc_to_tor(tor_grids::GridsT, ϕ_isl::Array{ComplexF64, 3}, isl_grids::GridsT, isl::IslandT)
    #TODO

    rgrid, θgrid, ζgrid = inst_grids(tor_grids)

    #may need to add the periodicity into this.
    κgrid, ᾱgrid, φgrid = inst_grids(isl_grids)

    #TODO
    #despite what Interpolations.jl says, we need to have the functions containing the periodic element
    #so fkn cooked
    #may just be worth getting the Hermite stuff working!
    isl_itp = interpolate((κgrid, ᾱgrid, φgrid), ϕ_isl, (Gridded(Linear()), Gridded(Linear(Periodic())), Gridded(Linear(Periodic()))));

    #this will obvs be a disaster for v values far from the island!
    isl_ext = extrapolate(isl_itp, Periodic());

    ϕ_tor = Array{ComplexF64}(undef, length(rgrid), length(θgrid), length(ζgrid));
    for (i, r) in enumerate(rgrid), (j, θ) in enumerate(θgrid), (k, ζ) in enumerate(ζgrid)

        κ, ᾱ, φ = tor_coords_to_isl(r, θ, ζ, isl)

        #originally, was 0, unsure exactly why this was here, or why it is needed.
        if κ==10000 && ᾱ==0 && φ==0
            #for some reason this case is being mapped to 2?
            ϕ_tor[i, j, k] = 0
        else

            ϕ_tor[i, j, k] = isl_ext(κ, ᾱ, φ)
        end
    end

    return ϕ_tor#, fft(ϕ_tor, [2, 3])

end

#anoying af that the arrays are abstract. This is because they are step ranges or whatever
function tor_efunc_to_isl!(ϕ_isl::Array{ComplexF64, 3}, κgrid::Array{Float64}, ᾱgrid::AbstractArray{Float64}, τgrid::AbstractArray{Float64}, ϕ_tor::Array{ComplexF64, 4}, rgrid::Array{Float64}, θgrid::AbstractArray{Float64}, ζgrid::AbstractArray{Float64}, isl::IslandT)

    for (i, κ) in enumerate(κgrid), (j, ᾱ) in enumerate(ᾱgrid), (k, τ) in enumerate(τgrid)

        #find the equivalent coordinate in toroidal coordinates.
        r, θ, ζ = isl_in_coords_to_tor(κ, ᾱ, τ, isl)

        #now we interpolate the value of ϕ_tor at this location.
        ϕ_isl[i, j, k] = hermite_interpolation(r, θ, ζ, ϕ_tor, rgrid, θgrid, ζgrid)

    end
end


#this version uses interpolation for the ϕ_tor.
function tor_efunc_to_isl!(ϕ_isl::Array{ComplexF64, 3}, κgrid::Array{Float64}, ᾱgrid::AbstractArray{Float64}, τgrid::AbstractArray{Float64}, ϕ_tor::Array{ComplexF64, 3}, rgrid::Array{Float64}, θgrid::AbstractArray{Float64}, ζgrid::AbstractArray{Float64}, isl::IslandT)

    #the periodic part of this is cooked.
    #can be fixed by periodifying the grids and ϕ, i.e. adding the first element to the end.
    tor_itp = interpolate((rgrid, θgrid, ζgrid), ϕ_tor, (Gridded(Linear()), Gridded(Linear(Periodic())), Gridded(Linear(Periodic()))));

    tor_ext = extrapolate(tor_itp, Periodic());

    #iterate through each point on our island grid
    for (i, κ) in enumerate(κgrid), (j, ᾱ) in enumerate(ᾱgrid), (k, τ) in enumerate(τgrid)

        #find the equivalent coordinate in toroidal coordinates.
        #fkn worst name so far.
        r, θ, ζ = isl_in_coords_to_tor(κ, ᾱ, τ, isl)

        #now we interpolate this coordinate using our functions.
        ϕ_isl[i, j, k] = tor_ext(r, θ, ζ)
        

    end
    
end

