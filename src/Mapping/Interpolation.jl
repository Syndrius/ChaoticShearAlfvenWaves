

function hermite_interpolation(r::Float64, θ::Float64, ζ::Float64, ϕ::Array{ComplexF64}, grids::FFFGridsT)

    #currently restricted to FFF
    #may also want a restriction of ϕ to ensure the derivative
    #although maybe that should be outside this function.

    #fkn stupid af to instantiate this each time....
    rgrid, θgrid, ζgrid = instantiate_grids(grids)
    #display(θgrid) 
    #display(ζgrid) 
     #assumes a uniform grid!!!
    #not true for r for fek sake.


    #these two are always uniform so we can just do this.
    dθ = θgrid[2] - θgrid[1]
    dζ = ζgrid[2] - ζgrid[1]

    #first we find the closest actual grid point.

    #so we actually need the two closest grid points...
    #think we can do this 1d at a time right???
    # I guess we will just ignore the derivative parts for now???
    # which is actually an awful idea I think
    #but we will have to change all our code to retrieve the derivatives....
    #we defs should be using the derivative...
    rind = find_ind(rgrid, r)
    θind = find_ind(θgrid, θ)
    ζind = find_ind(ζgrid, ζ)
    #ζind = 1

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
        #return ϕ[rind, θind, ζind]
        ##TODO
        if rind == 1
            r1 = rind
            r2 = rind+1
            dr = rgrid[r2] - rgrid[r1]
        elseif rind == grids.r.N
            r1 = rind
            #this is inconsistent with other cases, 
            #but should evaluate to zero for the other hermite function.
            r2 = rind-1
            dr = rgrid[r1] - rgrid[r2]
        else
            r1 = rind
            r2 = rind+1
            dr = rgrid[r2] - rgrid[r1]
        end
        #r2 = rind
        #r1 = rind
    elseif rgrid[rind] > r
        r1 = rind - 1
        r2 = rind
        dr = rgrid[r2] - rgrid[r1]
    else
        r1 = rind
        r2 = rind+1
        dr = rgrid[r2] - rgrid[r1]
    end

    #need some kind of periodicity here!
    if θgrid[θind] == θ 

        if θind == grids.θ.N #minus 1??
            θ1 = θind

            #inconsistent case but will evaluate to zero.
            θ2 = θind-1
        else
            θ1 = θind
            θ2 = θind + 1
        end
        #return ϕ[rind, θind, ζind, 1]
    elseif θgrid[θind] > θ
        θ1 = θind - 1
        θ2 = θind
    elseif θ > θgrid[end]
        #into the periodic zone
        θ1 = θind
        θ2 = 1
    else
        θ1 = θind 
        θ2 = θind + 1
    end


    if ζgrid[ζind] == ζ 

        if ζind == grids.ζ.N #minus 1??
            ζ1 = ζind

            #inconsistent case but will evaluate to zero.
            ζ2 = ζind-1
        else
            ζ1 = ζind
            ζ2 = ζind + 1
        end
        #return ϕ[rind, θind, ζind, 1]
    elseif ζgrid[ζind] > ζ
        ζ1 = ζind - 1
        ζ2 = ζind
    elseif ζ > ζgrid[end]
        #into the periodic zone
        ζ1 = ζind
        ζ2 = 1
    else
        ζ1 = ζind 
        ζ2 = ζind + 1
    end

    #TODO
    #need to figure out best way to do edge cases,
    #need to add derivs, need to add ζ
    #actually test this properly.
    

    #val2 = ϕ[r2, θind, ζind]
    #val1 = ϕ[r1, θind, ζind]

    Δr = (r - rgrid[r1]) / dr

    Δθ = (θ - θgrid[θ1]) / dθ

    Δζ = (ζ - ζgrid[ζ1]) /dζ

    #think we need to loop over ind pairs and each deriv.
    #perhaps we need to create an S.

    #not sure about the order of this...
    #S = rθζ, r'θζ, rθ'ζ, rθζ', r'θ'ζ, r'θζ', rθ'ζ', r'θ'ζ'
    #this should obviously not be defined here lol.
    grid_id = [0, 0, 1, 1]
    basis_id = [0, 1, 0, 1]


    rinds = [r1, r2]
    θinds = [θ1, θ2]
    ζinds = [ζ1, ζ2]
    ϕ_int = 0
    for hr in 1:4, hθ in 1:4, hζ in 1:4

        gi = (rinds[grid_id[hr]+1], θinds[grid_id[hθ]+1], ζinds[grid_id[hζ]+1])

        bi = 1 + basis_id[hζ] + 2 * basis_id[hθ] + 4 * basis_id[hr]
        #display(gi)
        ϕ_int += ϕ[gi..., bi] * hb(Δr, hr) * hb(Δθ, hθ) * hb(Δζ, hζ)
    end



    return ϕ_int
    

    #return ϕ[rind, θind, ζind]
end


#this should be inside basis
#and should be used lol
#not just arbitrary function definitions like we currently have.
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

#this should almost certainly be inside basis
function hb(x, ind) #think this needs to be kept in 1d actually!

    #big yucky.
    if ind == 1
        return h00(x)
    elseif ind == 2
        return h10(x)
    elseif ind == 3
        return h01(x)
    else
        return h11(x)
    end

end