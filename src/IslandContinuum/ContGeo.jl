
#mtric specifically for computing the island continuum.
mutable struct IslContMetT
    gl :: Array{Float64} #probably needs more info, no idea what size this will be!
    gu :: Array{Float64}
    J :: Array{Float64}
end


function isl_cont_metric(met, ψ, θ, geo) #what are these things!
    #we should be able to use our og function for this.
    #only diff is that that does things one grid point at a time
    #this does not, ideally we change this to match structure of other code.

    r = sqrt.(2*ψ ./geo.B0)

    dψdr = r .* geo.B0 #correct


    Δprime = r ./ (4*geo.R0)
    dΔprime = 1 / (4*geo.R0)

    #display(θ)
    ct = cos.(θ)
    st = sin.(θ)

    #display(ct)

    @. met.J = r * geo.R0 * (1 + 2 * r / geo.R0 * ct) / dψdr
    #display(metric.J)

    #these are actually the upper fellas, contrast to earlier...
    #not clear on the transformation from r to psi here...
    g¹¹ = @. dψdr^2 * (1+2*Δprime * ct) #so at least this one is wrong.

    #display(g¹¹)

    g¹² = @. - dψdr / r * st * (r/geo.R0 + Δprime + r*dΔprime)

    g²² = @. 1 / r^2 * (1 - 2 * (r/geo.R0 + Δprime) * ct)

    g³³ = @. 1 / geo.R0^2 /(1 + 2 * r / geo.R0 * ct)

    met.gu[1, 1, :, :] .= g¹¹
    #not sure if I will actually need both of these...
    met.gu[1, 2, :, :] .= g¹²
    met.gu[2, 1, :, :] .= g¹²
    met.gu[2, 2, :, :] .= g²²
    met.gu[3, 3, :, :] .= g³³

    #this could be the same as the jacobian....?
    scale_fact = @. -g¹²^2 * g³³ + g¹¹ * g²² * g³³

    #inverse is done manually...
    #this is awful!
    @. met.gl[1, 1, :, :] = g²² * g³³ / scale_fact
    #not sure if I will actually need both of these...
    @. met.gl[1, 2, :, :] = -g¹² * g³³ / scale_fact
    @. met.gl[2, 1, :, :] .= -g¹² * g³³ / scale_fact
    @. met.gl[2, 2, :, :] .= g¹¹ * g³³ / scale_fact
    @. met.gl[3, 3, :, :] .= (-g¹²^2 + g¹¹ * g²²) / scale_fact



end
