
#just going to compute the numerical integral for the mode overlap with a tae and possible (m_I, n_I) modes in the island.

#this file shows that there is at least some overlap between tae and generic island modes in island_damping_q case
#unclear which island modes actually overlap in frequency though!
#also hard to know if this is correct tbh.

using Elliptic
using Trapz

ψ0 = 0.125 #ie r=0.5
#start with ψ = 0.125, this will make κ's domain a bit simpler.
ψ = 0.125

#tae in form
#ϕ = ϕ_22(ψ)exp(i2θ - i2ζ) + ϕ_32(ψ)exp(i3θ - i2ζ)
#we don't really care about exponential term in the front, as tae is global, so there will be some radial overlap.

θgrid = range(0, 2π, 50)
ζgrid = range(0, 2π, 50)



#we probably need to know which island modes are above the og tae gap.

#begin with our original test case, ie 5/4 island.
m0 = 5
n0 = 4
A = 1e-4 #as a start 

qp = 0.4
q0 = m0/n0

#considering a single n in island contniuum case leads to no gap
#not sure how to split up different n's
m = 3
n = 2
#we can assume that ϕ_22 and ϕ_32 are the same, and set them to 1.

ϕ_tae = zeros(ComplexF64, length(θgrid), length(ζgrid))
#ie test of normal continuum case.
ϕ_cont = zeros(ComplexF64, length(θgrid), length(ζgrid))


#I guess a single island mode at a time?
ϕ_isl = zeros(ComplexF64, length(θgrid), length(ζgrid))

for (i, θ) in enumerate(θgrid), (j, ζ) in enumerate(ζgrid)

    ϕ_tae[i, j] = exp(2*1im*θ - 2 * 1im* ζ) + exp(3*1im*θ - 2 * 1im* ζ)

    ϕ_cont[i, j] = exp(2*1im*θ - 2 * 1im* ζ)
    #start by assuming ψ = ψ0
    κ = (1-cos(m0 * θ - n0 * ζ)) / 2
    width = 4 * sqrt(A * q0^2/qp)
    ψ̄ = 2*width/(m0*π) * ( (κ-1) * Elliptic.K(κ) + Elliptic.E(κ))

    ϕ_isl_coeff = exp(-ψ̄^2) # assume gaussian.

    asin_arg = 1/sqrt(κ) * sin((m0*θ-n0*ζ)/2)

    if abs(asin_arg) > 1 || κ==0
        #ᾱ = 0 #not sure about this... -> maybe a continue would be better?
        continue
    else
    #may need some if condition here tbh!

        ᾱ = π/(2*Elliptic.K(κ)) * Elliptic.F(asin(asin_arg), κ)
    end

    ϕ_isl[i, j] = ϕ_isl_coeff * exp(1im * m * ᾱ - 1im * n/m0 * ζ)


end


#ϕ_isl

#normal overlap of sin(x) and sin(x) is π, may need to compare with 
#regular continuum case?
#comparing overlap with cont case, typically 1-2 orders of magnitude smaller! may explain some things.
#we probably also need to compare the amplitudes of the island modes for more verification.
I = trapz((θgrid, ζgrid), ϕ_tae .* conj(ϕ_isl))
I = trapz((θgrid, ζgrid), ϕ_tae .* conj(ϕ_cont))
#I = trapz((θgrid, ζgrid), ϕ_res)

#res = 0.0 + 0.0im
#for (i, θ) in enumerate(θgrid), (j, ζ) in enumerate(ζgrid)

#    res += ϕ_isl[i, j] * ϕ_tae[i, j]
#end

#display(res)


