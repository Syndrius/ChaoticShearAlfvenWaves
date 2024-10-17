

function coords_tor_to_isl(r, θ, ζ, isl::IslandT)

    α = θ - ζ / isl.q0

    κ = 1/isl.w^2 * (r^2 - isl.r0^2)^2 + sin(isl.m0 * α / 2)^2

    if κ < 1
        #so this is at least partially wrong...
        τ = asin(1/sqrt(κ) * abs(sin(isl.m0 * α / 2)))

        if r > isl.r0
            #i.e. top right qudrant.
            #exact expressions need more explaining.
            #taken from restricted_mapping.jl, which uses abs to certify the correct quadrant,
            #may need to do that here as the 1/2 in the sin will probbaly cook things...
            if mod(α * isl.m0, 2π) < π

                ᾱ = π/ (2 * Elliptic.K(κ)) * Elliptic.F(τ, κ)

            else
                #top left quadrant, 4th by our clockwise notation.
                ᾱ = π/ (2 * Elliptic.K(κ)) * Elliptic.F(-τ + 2π, κ)
            end
        else

            if mod(α * isl.m0, 2π) < π
                #bottom right quadrant.
                ᾱ = π/ (2 * Elliptic.K(κ)) * Elliptic.F(-τ + π, κ)

            else
                #bottom left quadrant, 
                ᾱ = π/ (2 * Elliptic.K(κ)) * Elliptic.F(τ + π, κ)
            end
        end


        #ᾱ = π/ (2 * Elliptic.K(κ)) * Elliptic.F(τ, κ)
    else
        #this is still a work in progress.
        ᾱ = 0
        #without this we get some very surprising results, in that they match our real results pretty fkn well.
        return 0, 0, 0
    end

    return κ, ᾱ, ζ

end




function κ(r, θ, ζ, isl)

    

    #taking mod here cooks this apparently...
    #maybe a hint of what is going on??? Perhaps the elliptic stuff needs (-pi, pi)?
    #unsure why that would be????
    #α = mod(θ - ζ/isl.q0, 2π) #sign of this is questionable.
    α = θ - ζ/isl.q0 #sign of this is questionable.
    #χ = -isl.qp/(2*isl.q0^2) * (r^2/2 - r0^2/2)^2 + isl.A * cos(isl.m0*α)
    #this genuanly makes zero sense why this works, 
    #somewhere in this transformation or in our code we must be using r as psi...
    #χ = -isl.qp/(2*isl.q0^2) * (r - isl.r0)^2 + isl.A * cos(isl.m0*α)
    #χ = -isl.qp/(2*isl.q0^2) * (ψ - isl.r0)^2 + isl.A * cos(isl.m0*α)
    #κ = sqrt((-χ + isl.A) / (2*isl.A))
    #κ = ((-χ + isl.A) / (2*isl.A))
    κ = 1/isl.w^2 * (r^2 - isl.r0^2)^2 + sin(isl.m0 * α / 2)^2

    #display(κ)

    return κ
end


function κ_tor(κlist)

    #will ignore zeta for now
    rvals = zeros(length(κlist))
    θvals = zeros(length(κlist))
    #κvals = zeros(length(κlist), length(κlist))
    for κ in κlist
        α = 2/isl.m0 * asin(sqrt(κ)*Elliptic.Jacobi.sn(2*Elliptic.K(κ) * ᾱ / π, κ))

    end

    return rvals, θvals#, ζvals


end

        


