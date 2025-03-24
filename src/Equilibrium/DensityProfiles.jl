
#TODO change to flux!

function gae_isl_dens(r::Float64)
    κ = -1 #paper just says less than 0?
    p = 1
    fact = 5 #this is completly arbitrary, but just dicates the height of the peak.

    return fact*((1-p*r^2)^κ)
end


"""
    axel_dens(r::Float64)xs

Density function from Axel's paper, should be renamed.
"""
function axel_dens(r::Float64)
    #Axel's density function
    return 1/2 * (1 - tanh((r-0.8)/0.1))

end

"""
    uniform_dens(r::Float64)

Uniformd density function, default density, returns 1.0.
"""
function uniform_dens(r::Float64)
    return 1.0
end

"""
    bowden_singular_dens(r::Float64)

Density function from Bowden Singular paper, should be renamed.
"""
function bowden_singular_dens(r::Float64)
    return 1/2 * (1-tanh((r-0.7)/0.05))
end

"""
    comparison_bowden_dens(r::Float64)

Density function from Bowden comparison paper, should be renamed or removed as that paper seems dodge, only 800 points?.
"""
function comparison_bowden_dens(r::Float64)

    return (1-(r^2)^(7.5))^7
end
