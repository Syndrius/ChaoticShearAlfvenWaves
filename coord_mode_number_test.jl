
#testing mode numbers in different coords
using FFTW
using Plots; plotlyjs()
using Elliptic


x = range(0, 2π, 1000)

f(x) = real(exp(1im*x) + exp(2*1im*x))

plot(x, f.(x))


x_to_u(x) = x^2/(2*π)

u_to_x(u) = sqrt(2 * π * u)

#g(u) = real(exp(1im*u^2/(2*π)) + exp(2*1im*u^2/(2*π)))

g(x) = f(x_to_u(x))

u = sqrt.(2 .*π.*x)

plot(u, g.(u))

plot(u, f.(u) )

display(u .-x)

a = fft(g.(x))

z = zeros(ComplexF64, length(x))

for i in 1:800
    #z += @. a[i] * exp(-i * 1im * u_to_x(u))
    z += @. a[i] * exp(-i * 1im * x)
end

display(a[1:5])
#norm = maximum(abs.(z))

plot(u, real.((z)) )#, ylimits=(-5, 5))


#maybe try this with the actual transformation. bit unclear what is happening though!

m0 = 5
n0 = 4

A = 4e-4

θgrid = range(0, 2π, 50)
ζgrid = range(0, 2π, 50) #may not actually go this far!



#we will need to write this in terms of θ and ζ, not α.
αbar = zeros(length(θgrid), length(ζgrid))

#going to set ψ = ψ0

for i in 1:1:length(θgrid), j in 1:1:length(ζgrid)
    ang = m0 * θgrid[i] - n0 * ζgrid[j] #not sure about negative here.
    χ = cos(ang)
    #display(χ)
    κ = (-χ + 1)/(2)
    if κ < 1 #wot, this should be for inside the island??
        continue
    end
    αbar[i, j] = π / 2(Elliptic.K(κ)) * Elliptic.F(asin(1/sqrt(κ) * sin(ang/2)), κ)
end