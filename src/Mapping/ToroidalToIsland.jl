
#handles the requirement to map between from toriodal coordinates
#into island coordinates

function coords_isl_to_tor_Axel(κ, βs, φ, isl::IslandT)

    #w = 0.1 #fk me what the hel is width....

    #guess we can just try width the same as Zhisong? Probable has a factor of a half though!
    #and surely ι' is not equal to q'???
    #but if it is just a measure of the island width it should be ok right?
    
    #Axel uses half width, so we will take half the width of Zhisong.
    w = 2 * sqrt(isl.A * isl.q0^2 / isl.qp)

    #changing this just changes the distribution of the island, actually quite a minor effect tbh.
    #not sure what it correct.
    #κ = κ^2
    #κ = sqrt(κ)
    #I guess this cannot handle outside the ol island yet.
    if κ > 1 
        return 0, 0, 0
    end

    K = Elliptic.K(κ) 
    sinβ = Elliptic.Jacobi.sn(4*K * βs / 2π, κ)
    cosβ = Elliptic.Jacobi.cn(4*K * βs / 2π, κ)

    r = sqrt(isl.ψ0*2  + w*κ * cosβ)

    α = asin(κ * sinβ)

    θ = (2 * α - isl.n0 * φ)/isl.m0

    return r, θ, φ
end



#needs a different name by golly.
#we have had less success with this than others, may be worth ignoring for now.
function coords_isl_to_tor(κ, ᾱ, φ, isl::IslandT, sign)
    #name of this is a bit confusing
    #but idea is that we pass in island coordinates
    #and we find equivalent toroidal coordinates.

    if κ > 1
        α = 2/isl.m0 * Elliptic.Jacobi.am(isl.m0 * Elliptic.K(1/κ) * ᾱ / π, 1/κ)
        #return 0, 0, 0
    else
        α = 2/isl.m0 * asin(sqrt(κ)*Elliptic.Jacobi.sn(2*Elliptic.K(κ) * ᾱ / π, κ))
    end

    χ = -1 * (2*isl.A * κ - isl.A)

    #we need to be a careufl with how we have redefined χ, in particular there probbaly some factors 
    #of 1/2 going around from the r->ψ conversion.

    #may be a good idea to start writing the paper so that we have all our definitions written down tbh.

    #κ = 4/w^2(r^2/2 - r0^2/2)^2 + sin^@(m0 α / 2)
    # => κ = 1/w^2(r^2 - r0^2)^2 + sin^@(m0 α / 2)

    #lol wot even is r0...
    #taking abs is bold here, but only negatives are like e-20
    #res = sqrt(abs(-2*isl.q0^2/isl.qp * (χ - isl.A*cos(isl.m0 * α))))

    res = sqrt(isl.w^2 * (κ - sin(isl.m0*α/2)^2))

    αmod = mod(ᾱ, 2π) #probably not needed as we are passing in (0, 2π).

    θ = α + φ/isl.q0
    #for outside cases we have to choose a sign!
    if κ > 1
        if sign > 0
            r = sqrt(+res + isl.r0^2)
        else
            r = sqrt(-res + isl.r0^2)
        end

        #this doesn't seem to do anything??????
        #θ = mod(isl.m0 * θ, 2π)
        θ = mod(θ, 2π)
    else
        #inside the island both signs of ψ are needed.
        #r0 = sqrt(2*isl.ψ0)
        #not sure if ψ0 and ψ should be converted together or individiually.
        #think we probably need to do this differently now that the flux surfaces should be less cooked.
        if αmod < π/2
            r = sqrt(+res + isl.r0^2)
            #r = sqrt(2*res) + r0
        elseif αmod < 3π/2
            r = sqrt(-res + isl.r0^2)
            #r = -sqrt(2*res) + r0
        else
            r = sqrt(+res + isl.r0^2)
            #r = sqrt(2*res) + r0
        end
        θ = mod(θ, 2π)
    end

    
    #no fkn idea if this is actually flux or not lol.
    #display(res)
    #display((κ,  αmod, α))
    #if κ == 0.0
    #    r = 0.5
    #else
        #r = sqrt(2π)
    #end
    #return sqrt(2ψ), mod(θ, 2π), mod(φ, 2π) #may need an m0 here.
    #return ψ, mod(θ, 2π), mod(φ, 2π) #may need an m0 here.
    #return r, mod(θ, 2π), mod(φ, 2π) #may need an m0 here.
    return r, θ, mod(φ, 2π) #may need an m0 here.
    #return ψ + (0.5-0.125), mod(θ, 2π), mod(φ, 2π) #may need an m0 here.
    #return ψ, mod(θ, 2π), mod(φ, 2π)
    #return r, θ, φ
    #changing this along wtih ψ0 seems to have very minimal effect.
    #return ψ, θ, φ

end


#computes the sepratrix in toroidal coodinates.
function compute_sepratrix(grids, isl)

    sep1 = zeros(grids.θ.N)
    sep2 = zeros(grids.θ.N)

    #idealy we can generalise this.
    #qp = isl.qp
    #q0 = isl.q0
    #r0 = isl.r0

    _, θgrid, _ = inst_grids(grids)

    #κ = 4/w^2(r^2/2 - r0^2/2)^2 + sin^@(m0 α / 2)
    # => κ = 1/w^2(r^2 - r0^2)^2 + sin^@(m0 α / 2)

    #this is going to assume ζ=0 for simplicity, should generalise though!
    for i in 1:grids.θ.N
        α = θgrid[i]
        #unsure why this needs to be isl.A + isl.A -> should be negative???
        #res = sqrt(abs(-2 * q0^2 / qp * (isl.A + isl.A * cos(isl.m0 * α))))

        #should define this in terms of width rather than all of this garbage.
        #res = sqrt(abs(16 * q0^2 * r0 * isl.A / qp * (1 -  sin(isl.m0 * α/2)^2)))

        #sep is defined for κ=1
        res = sqrt(isl.w^2 * (1 - sin(isl.m0*α/2)^2))

        #ψsep1[i] = sqrt(2*(-res + 0.125))
        #ψsep2[i] = sqrt(2*(res + 0.125))
        #ψsep1[i] = -sqrt(2*(res)) + 0.5# + 0.125))
        #ψsep2[i] = sqrt(2*(res)) + 0.5# + 0.125))
        sep1[i] = sqrt(-res + isl.r0^2)
        sep2[i] = sqrt(res + isl.r0^2)
    end

    return sep1, sep2
end