





@kwdef struct ContIslandT
    m0 :: Int64 
    n0 :: Int64
    A :: Float64
    q0 :: Float64
    qp :: Float64
    ψ0 :: Float64 #ideally change this to r if possible
end

#this is cooked!
#don't want to be using these.
mutable struct B_field{A<:Array{Float64}, B<:Array{Float64}}
    B :: A
    B2 :: B
end



#just used for plotting
#p is for passing, probably need to be clearer with p for passing vs p for plus.
function ψ̄_p(isl, χ, sign)
    κ = @. (isl.A - χ) / (2*isl.A)
    w = 4 * sqrt(isl.A  * isl.q0^2 /isl.qprime)

    #may need a real for the elliptic fellas
    Ω = @. -sign * sqrt(isl.A * isl.qprime / isl.q0^2) * π * sqrt(κ) / real(ellipticK(1/κ))
    
    return @. isl.ψ0 + sign * w / π * sqrt(κ) * real(ellipticE(1/κ))
end


#start with the trapped particles

#looks to be working!
#may need to be able to handle lists of chi later though
function ω_t(χ, params)
    ω0 = sqrt(params.A * params.qprime / params.q0^2)

    κ = (-χ + params.A)/(2*params.A)

    K = real.(ellipticK(κ))
    return -ω0 * params.m0*π / (2*K) #does julia need returns?
end


#think the way this is done is terrible, involves much recomputation
#should just have a big function that computes everything, calling smaller functions as needed.
#need a better name for params here.
#matches python!
function α_t(χ, ᾱ, params)
    κ = (-χ + params.A)/(2*params.A)
    #println("kap")

    #display(K(κ))
    K = real.(ellipticK(κ))

    sn = @. real.(sin.(am.(2*K ./π .* ᾱ, κ))) #this is actually returing something wildly wrong for the last iteration!


    return @. 2/params.m0 * asin.(sqrt.(κ)*sn)

end

function ψ_t(χ, α, ᾱ, isl)

    #this is being recomputed in the same call of gradpsibar
    #α = α_t(χ, ᾱ, isl)


    #note division sign is floor division, equiv to // in python
    #not sure exactly why this form, but this tells which side of the island chain we are on, and hence the isgn of psi
    #may be able to use % somewhere in here to simplify this expr
    sign = 2*mod.( floor.((ᾱ .- π/2) / π), 2) .- 1

    #println(sign)
    #abs form is copied from python, not clear why it is in this form
    return @. sign * sqrt(2*isl.q0^2/isl.qprime) .* sqrt.(abs.( (isl.A * cos.(isl.m0*α) .- χ))) .+ isl.ψ0 
end

#don't think this is better for performance, but will be clearer what we are trying to do.
#think B is not a great name, as that is the struct that contains a B
function compute_Bfield!(B::B_field, ψ, α, isl, metric)


    #what the fuk theta grid is passed in here? wot is going on!!!!
    #this is a very stupid function,
    #we need to re-broadcast theta and zeta so that this will have the same shape as the metric, which is computed
    #based on psi, which is computed from, and has the shape of, alphabar, seems like madness to do it this way.
    #feels like we should just craft some kind of coordinate array or something, instead of recreating it all over the place.

    #I think this is the pre new coords magnetic field, really not clear though!

    #needs to be broadcast I think!
    q = @. 1 / (1 / isl.q0 - isl.qprime /isl.q0^2 * (ψ - isl.ψ0))

    #Bζ = 1 ./ metric.J
    #Bθ = 1 ./ (metric.J .* q)
    #Bψ = @. isl.A * isl.m0 *sin.(isl.m0 .* α)

    #display(Bψ)
    #this needs to have been done earlier!
    #B = zeros(3, size(ψ)...)
    #this may need to be tranformed in some way for this to work as expected.
    @. B.B[1, :, :] = isl.A * isl.m0 *sin(isl.m0 * α)
    @. B.B[2, :, :] = 1 / (metric.J * q)
    @. B.B[3, :, :] = 1 / metric.J 


    #also should have been done earlier!
    B2 = zeros(size(ψ)...) #not sure this stuff is needed, because this code always assumes a 2d grid.


    #not ideal!
    #not sure if this is a bottleneck, but this really needs to be fixed.
    #cannot be confident this works in general, but this does seem to work for test case at least.
    for i in 1:1:size(ψ)[1]
        for j in 1:1:size(ψ)[2]
            for k in 1:1:3

                for l in 1:1:3

                    B2[i, j] += B.B[k, i, j] * metric.gl[k, l, i, j] * B.B[l, i, j]
                end
            end
        end
    end
    #clunky solution, should be a better way to do the above matrix multiplication!
    B.B2 .= B2
    #println("B2")
    #println(B2)
    #println("Stop")

    #return B, B2

end

function ∇ψ̄2_t!(∇ψ̄2, χ, ψ, α, metric, isl)

    ω2 = ω_t(χ, isl)^2

    @. ∇ψ̄2 =  1/ω2 * (isl.qprime ^2 /isl.q0^4 * (ψ - isl.ψ0)^2 * metric.gu[1, 1, :, :] +
                isl.A^2 * isl.m0 *sin(α *isl.m0)^2 * (metric.gu[2, 2, :, :] + metric.gu[3, 3, :, :] / isl.q0^2) +
                2*isl.qprime / isl.q0^2 * (ψ - isl.ψ0) * isl.A * isl.m0 *sin(isl.m0 *α) * metric.gu[1, 2, :, :])
end


function ω_p(χ, params, sign)
    ω0 = sqrt(params.A * params.qprime / params.q0^2)

    κ = (-χ + params.A)/(2*params.A)

    
    return -sign*ω0 * π * sqrt(κ)/ellipticK(1/κ)
end


#think the way this is done is terrible, involves much recomputation
#should just have a big function that computes everything, calling smaller functions as needed.
#need a better name for params here.
#matches python!
function α_p(χ, ᾱ, params)
    κ = (-χ + params.A)/(2*params.A)
    #println("kap")
    print(κ)
    #κ should always be bigger than 1 outside the island, 
    #so 1/κ should be fine, I guess for the no island case we are getting some cookery?
    K = real.(Elliptic.K(1/κ))

    #display(params.m0*K .* ᾱ ./π)
    #display(1/κ)

    #whatever this am function is seems to be cooked for just the am, fine when you take the sin tho??
    #amval = am.(params.m0*K .* ᾱ ./π, 1/κ)
    #this form of am needs a real number though!
    amval = @. Jacobi.am.(params.m0*K * ᾱ /π, 1/κ)
    #display(amval)

    return @. 2 / params.m0 * amval

end

function ψ_p(χ, α, isl, sign)

    #this is being recomputed in the same call of gradpsibar
    #α = α_p(χ, ᾱ, isl)

    #display(α)

    #println(sign)
    #abs form is copied from python, not clear why it is in this form
    return @. sign * sqrt(2*isl.q0^2/isl.qprime) .* sqrt.(isl.A * cos.(isl.m0*α) .- χ) .+ isl.ψ0 
end

#not sure if this really needs to be differenent in the t or p case.
function ∇ψ̄2_p!(∇ψ̄2, χ, ψ, α, metric, isl, sign)

    ω2 = ω_p(χ, isl, sign)^2 #looks like we don't need this anymore...
    ω2 = 1 #think we are now actually getting chi not psibar??

    #display(metric.gu[1, 1, :, :])
    #display(metric.J)

    @. ∇ψ̄2 =  1/ω2 * (isl.qprime ^2 /isl.q0^4 * (ψ - isl.ψ0)^2 * metric.gu[1, 1, :, :] +
                isl.A^2 * isl.m0 *sin(α *isl.m0)^2 * (metric.gu[2, 2, :, :] + metric.gu[3, 3, :, :] / isl.q0^2) +
                2*isl.qprime / isl.q0^2 * (ψ - isl.ψ0) * isl.A * isl.m0 *sin(isl.m0 *α) * metric.gu[1, 2, :, :])
end


function ψ_p_new(χ, α, isl, sign)

    #this is being recomputed in the same call of gradpsibar
    #α = α_p(χ, ᾱ, isl)

    #display(α)

    #println(sign)
    #abs form is copied from python, not clear why it is in this form
    return @. sign * sqrt(2*isl.q0^2/isl.qp) .* sqrt.(isl.A * cos.(isl.m0*α) .- χ) .+ isl.ψ0 
end
