

using MID
using Elliptic
using Plots; plotlyjs()
using FFTW

#mayhap we should start with n=0???

#so we need to map our grid into (ψ, α, ζ coords?) unsure what this will look like.
#or do we pick a new grid and map from before??/

#genuinly have no idea what the fk Axel's κ is supposed to be.
#I think we should just try with Zhisongs ψ̄ first.

#og grids.
Nr = 10
Nθ = 5
Nζ = 4
rgrid = init_fem_grid(N=Nr)
θgrid = init_fem_grid(N=Nθ)
ζgrid = init_fem_grid(N=Nζ)

grids = init_grids(Nr, θgrid, ζgrid)

#unsure how big to make the new grids??
#new coords
Nκ = 157 #use kappa instead, just skips some steps I think.
Nα = 34 #feel like this one especially is completly cooked 
Nφ = 21



#ok perhaps an easier approach is to define a simple wave inside an island in island coords, and map to normal coords to see what it looks like?

#so ψ̄ is a function of kappa so we should be able to map to kappa for plotting as comparison...
#ideally we want a sin/cos in ψ to reflect Axel's modes, combined with m=1 etc in α, think we can just ignore φ for now.
ϕ_isl = zeros(ComplexF64, Nκ, Nα, Nφ);

#\kappa going to 2 is kinda fake, but should work for our purposes I think.
#these should be barred!!!!
κgrid = LinRange(0, 2, Nκ)
αgrid = LinRange(0, 2π, Nα)#+1)[1:end-1]
φgrid = LinRange(0, 2π, Nφ)#+1)[1:end-1]

function test_mode(κ, α, φ)
    #so kappa is 0<κ<1
    #α between 0 and 2π
    #may need a fkn m0 in here somehwere!
    if κ > 1
        return 0
    end

    m=1

    #return exp(-(κ-0.9)^2/0.01) * exp(1im * m * α)

    return (-4*(κ-0.5)^2+1)*exp(1im * m * α)
    #return (-4*(κ-0.5)^2+1)*cos(m * α)

    #try regular continuum mode...

    if κ < 1
        return 0
    end

    return exp(-(κ-4.7)^2/0.01) * exp(1im * m * α)

end

for (i, κ) in enumerate(κgrid), (j, α) in enumerate(αgrid), (k, φ) in enumerate(φgrid)
    ϕ_isl[i, j, k] = test_mode(κ, α, φ)
end

surface(αgrid, κgrid, real.(ϕ_isl[:, :, 1]))

Nr = 191
Nθ = 67
Nζ = 23

#this did not fix the island lack of sym
#ψgrid = LinRange(0, 0.5, Nr)
#rgrid = sqrt.(2 .*ψgrid)

rgrid = LinRange(0, 1, Nr)

θgrid = LinRange(0, 2π, Nθ)#+1)[1:end-1]
ζgrid = LinRange(0, 2π, Nζ)#+1)[1:end-1]
#ok so can we map this to a normal ϕ??? 
ϕ_normy = zeros(ComplexF64, Nr, Nθ, Nζ);

function coord_map(r, θ, ζ)
    ψ = r^2/2

    qp = 1.6
    m0 = 4
    n0 = 2
    q0 = m0/n0
    A = 1e-3
    ψ0 = 0.125
    #minus here seems to make no difference...
    α = θ - ζ/q0 #this minus will probably cook us lol.
    χ = -qp/(2*q0^2) * (ψ - ψ0)^2 + A * cos(m0 * α)

    κ = (-χ + A) / (2 * A) #so I guess we just ignore ψ̄???

    if κ > 1 #outside the island
        #second m0 is probably a 2.
        #cannot tell this from this case as we are zero outside island.
        #ᾱ = Elliptic.F(m0/2 * α, 1/κ) * π/(m0 * Elliptic.K(1/κ))
        #ᾱ = rem(Elliptic.F(m0/2 * α, 1/κ) * π/(2 * Elliptic.K(1/κ)), 2*m0*π)
        #ᾱ = Elliptic.F(m0/2 * α, 1/κ) * π/(m0 * Elliptic.K(1/κ))
        ᾱ = Elliptic.F(m0/2 * α, 1/κ) * π/(m0 * Elliptic.K(1/κ)) #+ π/2


    else
        #this can error, is this a rounding error
        #or are there more specific values for α???
        #this seems like a problematic solution.
        #only cases seen so far give asin(1 + ϵ). so perhaps it is just some rounding garbage.
        arg = 1/sqrt(κ)*sin(m0*α/2)
        try
            ᾱ = π/(2*Elliptic.K(κ)) * Elliptic.F(asin(arg), κ)
        catch
            #
            #display((κ, α))
            #display((1/sqrt(κ)*sin(m0*α/2)))
            #display(κ)
            #display(arg)

            #this corresponds to the island centre's
            #first island seems to hit this trigger without this if condition, may just be numerical coincidence, but there is other problemos with that island as well.
            if κ == 0.0
                ᾱ=NaN
            
            elseif arg < 0
            #if arg < 0
                ᾱ = π/(2*Elliptic.K(κ)) * Elliptic.F(π/2, κ)
            else
                ᾱ = π/(2*Elliptic.K(κ)) * Elliptic.F(3*π/2, κ)
            end
        end
        #this makes a big difference, function is now a smooth circle
        #but it is also a smooth circle for m=2 island mode...
        #nope wasn't this..., 
        #had a div
        #ᾱ = rem(ᾱ, 2π)
        #ᾱ -= π
    end
    #if ᾱ==NaN
    #    display(r, θ, ζ)
    #end

    return κ, ᾱ, ζ

end

function coord_to_ind(κ, α, φ)

    #this will need to be more sophisticated when we consider Hermite splines
    #this should be adequate for now though.
    #think the mod fixes this.
    #this is because α can be negative.
    #still seems kind of sub optimal...
    #return argmin(abs.(κgrid .- κ)), argmin(abs.(αgrid .- mod(α, 2π))), argmin(abs.(φgrid .- φ))
    return argmin(abs.(κgrid .- κ)), argmin(abs.(αgrid .- α)), argmin(abs.(φgrid .- φ))
    #return argmin(abs.(κgrid .- κ)), argmin(abs.(αgrid .- abs(α))), argmin(abs.(φgrid .- φ))
end


albar_vals = zeros(Float64, Nr, Nθ, Nζ);

for (i, r) in enumerate(rgrid), (j, θ) in enumerate(θgrid), (k, ζ) in enumerate(ζgrid)

    #going to need more inputs.
    κ, α, φ = coord_map(r, θ, ζ)
    #if α < 0
    #    println("Hello")
    #end
    #this div here makes a huge difference!
    #makes it look good for m=1 case, but not for any other....
    #something is wrong with the output alphabar values...
    #κind, αind, φind = coord_to_ind(κ, div(α, π), φ)
    #think this mod is required, maps α to (0, 2π), otherwise it just maps to the zero point.
    κind, αind, φind = coord_to_ind(κ, mod(α, 2π), φ)
    ϕ_normy[i, j, k] = ϕ_isl[κind, αind, φind]
    albar_vals[i, j, k] = α

    #if isnan(α) && κ==0.0
    #    display((κ, α, φ))
    #    display((r, θ, ζ))
    #end

end
contourf(θgrid, rgrid, real.(albar_vals[:, :, 1]), levels=50)

#display(ζgrid)
#Still seems to be a small issue at θ=0, 2π, this may just be because the last island has its centre here though!
#so the surface seem to work when ᾱ = -π/2 but not for the positive π/2 side... wot. I guess this is because our arbitrary ᾱ grid is from 0 to 2π... 

#think we have fixed most issues now, but we still have the up-down symmetry problemo..
#inverse looks ok, still not v smooth
#mode outside the island is more cooked than anything lol
#especially case with A=1e-6, completly bananas.
surface(θgrid, rgrid, real.(ϕ_normy[:, :, 1]))

contourf(θgrid, rgrid, real.(ϕ_normy[:, :, 1]), levels=50)

display(θgrid[5])
cos(3 * (θgrid[5]-(2/3)*ζgrid[3]))
#think this is highlihgting the problemo.
#the values don't make any sense
#island should also be syymetric surely???
#also island split sorta works for 2 of them, not for the third, we are getting an extra nan in one of the centres, which is cooked af.
contourf(θgrid, rgrid, rem.(real.(albar_vals[:, :, 1]), π))
contourf(θgrid, rgrid, real.(albar_vals[:, :, 1]), levels=50)

minimum(albar_vals)
#ok, now can we turn this back into an island mode???
#then we can, in theory, do some actual modes.


Nκ2 = 194 #use kappa instead, just skips some steps I think.
Nα2 = 141 #feel like this one especially is completly cooked 
Nφ2 = 27



#ok perhaps an easier approach is to define a simple wave inside an island in island coords, and map to normal coords to see what it looks like?

#so ψ̄ is a function of kappa so we should be able to map to kappa for plotting as comparison...
#ideally we want a sin/cos in ψ to reflect Axel's modes, combined with m=1 etc in α, think we can just ignore φ for now.


#\kappa going to 2 is kinda fake, but should work for our purposes I think.
#these should be barred!!!!
κgrid2 = LinRange(0, 5, Nκ2)
αgrid2 = LinRange(0, 2π, Nα2)
φgrid2 = LinRange(0, 2π, Nφ2)

ϕ_isl2 = zeros(ComplexF64, Nκ2, Nα2, Nφ2);
function coord_map2(κ, ᾱ, φ)
    #a little bit unsure how to deal with inside vs outside.
    ζ = φ

    m0 = 3
    n0 = 2
    q0 = m0/n0
    qp = 1.6
    ψ0 = 0.125
    A = 1e-4

    χ = -(κ*2*A - A)

    

    if κ > 1
        #unsure about this...
        #this may be different outside.
        #hard to know though, but I think the fourier rep is different
        #outside vs inside, even for this case???
        #α = 2/m0 * Elliptic.Jacobi.am(m0*Elliptic.K(1/κ)/π, 1/κ)
        α = 2/m0 * Elliptic.Jacobi.am(2*Elliptic.K(1/κ)/π * ᾱ, 1/κ)
        θ = α + ζ/q0
        ψ = sqrt(-2*q0^2/qp * (χ - A * cos(m0 * α))) + ψ0
        #ψ = 0
    else
        α = 2/m0 * asin(sqrt(κ)*Elliptic.Jacobi.sn(2*Elliptic.K(κ)/π * ᾱ, κ))
        #this seems to make very little difference...
        #unsure what this should be....
        θ = α + ζ/q0
        arg = -2*q0^2/qp * (χ - A * cos(m0 * α))
        if arg < 0
            #this is triggering a lot, but always for e-20, so irrelevant
            #ideally we could tell julia to ignore this.
            #display(arg)
            arg = abs(arg)
        end
        ψ = sqrt(arg) + ψ0
        #ψ = 0.125
    end
    #this will obvs cause issues.
    

    r = sqrt(2*ψ)

    return r, θ, ζ
end

function coord_to_ind2(r, θ, ζ)

    return argmin(abs.(rgrid .- r)), argmin(abs.(θgrid .- mod(θ, 2π))), argmin(abs.(ζgrid .- mod(ζ, 2π)))
    #return argmin(abs.(rgrid .- r)), argmin(abs.(θgrid .- θ)), argmin(abs.(ζgrid .- ζ))
end


for (i, κ) in enumerate(κgrid2), (j, α) in enumerate(αgrid2), (k, φ) in enumerate(φgrid2)   

    r, θ, ζ = coord_map2(κ, α, φ)
    rind, θind, ζind = coord_to_ind2(r, θ, ζ)
    ϕ_isl2[i, j, k] = ϕ_normy[rind, θind, ζind]
end

#mayhap this is correct??? v hard to tell.
#could be a mistake, or could just be that this is quite numerically unstable
#increase of res makes this a lot better, still a bit off though, could be some minor errors somewhere
#and we are just getting the closest points so could be cooked af
#transforming to normal coords gives a pretty smooth function though???
#one issue could be the θ vs α grid, and the 2π vs 2m0π.
#think the weird spikes for κ > 1 are a real concern...
surface(αgrid2, κgrid2, real.(ϕ_isl2[:, :, 1]))

#still need to get it into normal form.
ϕ_fft2 = fft(ϕ_isl2, [2, 3]);
ϕ_fft = fft(ϕ_isl, [2, 3]);

p = plot(xlimits=(-0.05, 1.05))

for i in 1:Nα2
    plot!(κgrid2, real.(ϕ_fft2[:, i, 1]))
end
display(p)

p = plot(xlimits=(-0.05, 1.05))
for i in 1:Nα
    plot!(κgrid, real.(ϕ_fft[:, i, 1]))
end
display(p)


ϕ_fft_noi = fft(ϕ_normy, [2, 3]);

p = plot()

for i in 1:Nθ
    plot!(rgrid, real.(ϕ_fft_noi[:, i, 3]))
end
#this has an unbelivably clean ft, wot. How is that possible???
#like a perfect m=-1 ft, with the same double bump shape we have previously seen???? wild.
#seems a bit to good to be true, kinda sketch lol.
#would be nice to be able to prove this case, cause this is wild???
#how the fk does this happen?
#then non=n=0 cases are a bit more complicated.
#still kind of wild!!!
#the n=0 mode is an order of magnitude above the others, but perhaps not as dominant as we might expect.
#the n=2 mode is hude, i.e. on par with the n=0. seems like much of the complexity is inside the toroidal modes rather than poloidal???
#bit wild.
#same as the n=4, mode, i guess with n0=2, everysecond mode is quite dominant.
#interesting we don't see this behaviour in our code, at least not this prominantly.
display(p)
#seeing vague hints of og structure, 
#in partcular, the main mode is (one of) the most prominant.
#interestingly, the m=0 mode is doing weird stuff, especially at the island centre???
#we also note that the behaviour around the sepratrix has been modified
#to be expected I think!
#think this is a good start, we obviously need to do this
#for the actual case with Hermite bois.
#and I am a bit unclear about the size of the new grids, I guess they are just arbitrary??? I suppose we just keep increasing them until nothing changes???

#seams very odd that the transformation one way is smooth as butter, but then blocky and kinda garbage the other way, perhaps this is just how the coordinates happen to be spread??? -> or its because we are mapping a tiny r  in (0.4, 0.6) region into the full area, so obvs will be blocky af.
#hoepfully with our clustering this is not a problemo.
#may have to test this with higher resolution in this case.


#Qs on this. Transformation inside vs outside
#size of resulting grid
#use of κ
#m0 in appendix for outside
#ψ0 gone for inside
#taking mod of angles
#why no plus minus for 
#need to do this for a continuum mode!!!!!
#i.e. localised to a specific κ