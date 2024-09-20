
#handles the requirement to map between from toriodal coordinates
#into island coordinates

function coords_isl_to_tor(κ, βs, φ, isl::ContIslandT)

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
function coords_isl_to_tor_qu(κ, ᾱ, φ, isl::ContIslandT)
    #name of this is a bit confusing
    #but idea is that we pass in island coordinates
    #and we find equivalent toroidal coordinates.

    if κ > 1
        α = 2/isl.m0 * Elliptic.Jacobi.am(isl.m0 * Elliptic.K(1/κ) * ᾱ / π, 1/κ)
    else
        α = 2/isl.m0 * asin(sqrt(κ)*Elliptic.Jacobi.sn(2*Elliptic.K(κ) * ᾱ / π, κ))
    end

    χ = -1 * (2*isl.A * κ - isl.A)

    #lol wot even is r0...
    #taking abs is bold here, but only negatives are like e-20
    ψ = sqrt(abs(-2*isl.q0^2/isl.qp * (χ - isl.A*cos(isl.m0 * α)))) + isl.ψ0

    θ = α + φ/isl.q0
    #no fkn idea if this is actually flux or not lol.
    #return sqrt(2*ψ), θ, φ
    #changing this along wtih ψ0 seems to have very minimal effect.
    return ψ, θ, φ

end