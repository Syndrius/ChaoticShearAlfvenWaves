
#can we do some hermite interpolation???
#come on down...
using MID
using FFTW
using Plots; plotlyjs()

Nr = 30
Nθ = 6
Nζ = 1
rgrid = init_fem_grid(N=Nr);
θgrid = init_fem_grid(N=Nθ, pf=2);
ζgrid = init_fem_grid(N=Nζ, pf=-2);
θgrid = init_sm_grid(start=2, count = 2)
ζgrid = init_sm_grid(start=-2, count=1);

grids = init_grids(rgrid, θgrid, ζgrid);

geo = GeoParamsT(R0=10.0)

prob = init_problem(q=Axel_q, geo=geo); 


#with @views. 22.907080 seconds (8.07 M allocations: 721.379 MiB, 1.36% gc time)
#fk load more allocations and gc without views.
#outrageous spead up shifting the ϕ[:, test, :, ...] to ϕ[:, testr, testθ, :, :]

evals, ϕ, ϕft = compute_spectrum(prob=prob, grids=grids, full_spectrum=true)#, deriv=true); 


plot_continuum(evals)


ind = find_ind(evals, 0.383)
#ind = 348
plot_potential(ϕft, grids, ind)

plot_potential(ϕ, grids, ind)

contour_plot(ϕ, grids, ind=ind)

#test_fn = ϕ[ind, :, :, :, :];
test_fn = ϕ[ind, :, :, :];


function S(x, y, z) #think this needs to be kept in 1d actually!

end


function interpolate_phi(r, θ, ζ, ϕ, grids)

    rgrid, θgrid, ζgrid = instantiate_grids(grids)
    #display(θgrid) 
    #display(ζgrid) 
    dr = rgrid[2] - rgrid[1]
    dθ = θgrid[2] - θgrid[1]
    #dζ = ζgrid[2] - ζgrid[1]

    #first we find the closest actual grid point.

    #so we actually need the two closest grid points...
    #think we can do this 1d at a time right???
    # I guess we will just ignore the derivative parts for now???
    # which is actually an awful idea I think
    #but we will have to change all our code to retrieve the derivatives....
    #we defs should be using the derivative...
    rind = find_ind(rgrid, r)
    θind = find_ind(θgrid, θ)
    #ζind = find_ind(ζgrid, ζ)
    ζind = 1

    #may want some kind of check in case it is exact???
    #this works ok, but it looks like we need derivatives... RIP.
    #first lets see if we can do this into multiple dimensions.
    if rgrid[rind] == r

        #r1 = rind
        #r2 = 0
        #this only works because we are only interpolating in 1d
        #will be trickier for general case
        #do we just need like 7 if conditions??
        #surely there is a better way to do this??
        return ϕ[rind, θind, ζind]
        ##TODO
        if rind == 1
            r1 = rind
            r2 = rind+1
        elseif rind == grids.r.N
            r1 = rind
            #this is inconsistent with other cases, 
            #but should evaluate to zero for the other hermite function.
            r2 = rind-1
        else
            r1 = rind
            r2 = rind+1
        end
        #r2 = rind
        #r1 = rind
    elseif rgrid[rind] > r
        r1 = rind - 1
        r2 = rind
    else
        r1 = rind
        r2 = rind+1
    end

    #need some kind of periodicity here!
    if θgrid[θind] == θ || θind==6
        return ϕ[rind, θind, ζind, 1]
    elseif θgrid[θind] > θ
        θ1 = θind - 1
        θ2 = θind
    else
        θ1 = θind 
        θ2 = θind + 1
    end

    #TODO
    #need to figure out best way to do edge cases,
    #need to add derivs, need to add ζ
    #actually test this properly.
    

    #val2 = ϕ[r2, θind, ζind]
    #val1 = ϕ[r1, θind, ζind]

    Δr = (r - rgrid[r1]) / dr

    Δθ = (θ - θgrid[θ1]) / dθ

    #think we need to loop over ind pairs and each deriv.
    #perhaps we need to create an S.

    #not sure about the order of this...
    #S = rθζ, r'θζ, rθ'ζ, rθζ', r'θ'ζ, r'θζ', rθ'ζ', r'θ'ζ'

    #this is going to take too much thinking to do now!
    ϕ_int = 0
    for ri in [r1, r2], θi in [θ1, θ2], rd in [0, 1], θd in [0, 1]
        ϕ_int += 1

    end
    #think this is the way. come up with all combos of this
    #then put it into a loopy boi.
    ϕ_int += ϕ[r1, θ1, ζind, 1] * h00(Δr) * h00(Δθ)



    ϕ_int = ϕ[r1, θ1, ζind] * h00(Δr) * h00(Δθ) + ϕ[r2, θ1, ζind] * h01(Δr) * h00(Δθ) + ϕ[r1, θ2, ζind] * h00(Δr) * h01(Δθ) + ϕ[r2, θ2, ζind] * h01(Δr) * h01(Δθ)
    #think deriv should be 2. unsure though!
    #ϕ_int = ϕ[r1, θind, ζind, 1] * h00(Δr)  + ϕ[r2, θind, ζind, 1] * h01(Δr) #+ ϕ[r1, θind, ζind, 2] * h10(Δr)  + ϕ[r2, θind, ζind, 2] * h11(Δr)
    #ϕ_int = ϕ[r1, θind, ζind, 1] * h00(Δr)  + ϕ[r2, θind, ζind, 1] * h01(Δr) #+ ϕ[r1, θind, ζind, 2] * h10(Δr)  + ϕ[r2, θind, ζind, 2] * h11(Δr)
    #ϕ_int = ϕ[r1, θind, ζind] * h00(Δr)  + ϕ[r2, θind, ζind] * h01(Δr) #+ ϕ[r1, θind, ζind, 2] * h10(Δr)  + ϕ[r2, θind, ζind, 2] * h11(Δr)
    #* h00(Δθ) + ϕ[r1, θ2, ζind] * h00(Δr) * h01(Δθ) + ϕ[r2, θ2, ζind] * h01(Δr) * h01(Δθ)

#    display(Δr)


    return ϕ_int
    

    #return ϕ[rind, θind, ζind]
end


#can we interpolate this bad boi???
Nr = 100
Nθ = 30
Nζ = 15

rtest = LinRange(0, 1, Nr)
θtest = LinRange(0, 2π, Nθ+1)[1:end-1]
ζtest = LinRange(0, 2π, Nζ+1)[1:end-1]

ϕ_int = zeros(ComplexF64, Nr, Nθ, Nζ);
for (i, r) in enumerate(rtest), (j, θ) in enumerate(θtest), (k, ζ) in enumerate(ζtest)

    ϕ_int[i, j, k] = interpolate_phi(r, θ, ζ, test_fn, grids)


end


ϕ_int_fft = fft(ϕ_int, [2, 3]);

p = plot()
for i in 1:Nθ
    plot!(rtest, real.(ϕ_int_fft[:, i, 1]))
    #plot!(rtest, real.(ϕ_int[:, i, 1]))
end
display(p)

#This is not working lol.
#seems to be much worse than before...
#needs θ interp as well for this to look anything like the og tae
#can cleary see the fkery for θ > grids.θ.N
#think we are making progress but not heaps lol.
contourf(θtest, rtest, real.(ϕ_int[:, :, 1]), color=:turbo, levels=100)


function h00(t)

    return (1+2*t)*(1-t)^2
end

function h10(t)
    return t*(1-t)^2
end

function h01(t)
    return t^2*(3-2*t)
end

function h11(t)
    return t^2*(t-1)
end

t = LinRange(0, 1, 100)

#this does make sense, if each neighbouring grid points have the same value, the inbetween will also have the same value, unless the gradient is different.
#plot(t, h00.(t) .+ 2 .* h01.(t))
plot(t, h00.(t))
plot!(t, h01.(t))
plot!(t, h10.(t))
plot!(t, h11.(t))
