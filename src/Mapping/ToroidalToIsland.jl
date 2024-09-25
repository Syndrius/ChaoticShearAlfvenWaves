
#handles the requirement to map between from toriodal coordinates
#into island coordinates

function coords_isl_to_tor_Axel(κ, βs, φ, isl::ContIslandT)

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
function coords_isl_to_tor(κ, ᾱ, φ, isl::ContIslandT)
    #name of this is a bit confusing
    #but idea is that we pass in island coordinates
    #and we find equivalent toroidal coordinates.

    if κ > 1
        α = 2/isl.m0 * Elliptic.Jacobi.am(isl.m0 * Elliptic.K(1/κ) * ᾱ / π, 1/κ)
        return 0, 0, 0
    else
        α = 2/isl.m0 * asin(sqrt(κ)*Elliptic.Jacobi.sn(2*Elliptic.K(κ) * ᾱ / π, κ))
    end

    χ = -1 * (2*isl.A * κ - isl.A)



    #lol wot even is r0...
    #taking abs is bold here, but only negatives are like e-20
    res = sqrt(abs(-2*isl.q0^2/isl.qp * (χ - isl.A*cos(isl.m0 * α))))

    αmod = mod(ᾱ, 2π) #probably not needed as we are passing in (0, 2π).
    #r0 = sqrt(2*isl.ψ0)
    #not sure if ψ0 and ψ should be converted together or individiually.
    if αmod < π/2
        ψ = +res + isl.ψ0
        #r = sqrt(2*res) + r0
    elseif αmod < 3π/2
        ψ = -res + isl.ψ0
        #r = -sqrt(2*res) + r0
    else
        ψ = +res + isl.ψ0
        #r = sqrt(2*res) + r0
    end

    θ = α + φ/isl.q0
    #no fkn idea if this is actually flux or not lol.
    #display(res)
    #display((κ,  αmod, α))
    #if κ == 0.0
    #    r = 0.5
    #else
        #r = sqrt(2π)
    #end
    return sqrt(2ψ), mod(θ, 2π), mod(φ, 2π) #may need an m0 here.
    return ψ, mod(θ, 2π), mod(φ, 2π) #may need an m0 here.
    #return ψ + (0.5-0.125), mod(θ, 2π), mod(φ, 2π) #may need an m0 here.
    #return ψ, mod(θ, 2π), mod(φ, 2π)
    #return r, θ, φ
    #changing this along wtih ψ0 seems to have very minimal effect.
    #return ψ, θ, φ

end