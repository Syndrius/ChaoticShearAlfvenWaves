


"""
Density function from Axel's paper, should be renamed.

# Args
r::Float64 - Radial point.
"""
function axel_dens(r::Float64)
    #Axel's density function
    return 1/2 * (1 - tanh((r-0.8)/0.1))

end

"""
UniformdDensity function, default density.

# Args
r::Float64 - Radial point.
"""
function uniform_dens(r::Float64)
    return 1.0
end

"""
Density function from Bowden Singular paper, should be renamed.

# Args
r::Float64 - Radial point.
"""
function bowden_singular_dens(r::Float64)
    return 1/2 * (1-tanh((r-0.7)/0.05))
end

"""
Density function from Bowden comparison paper, should be renamed or removed as that paper seems dodge, only 800 points?.

# Args
r::Float64 - Radial point.
"""
function comparison_bowden_dens(r::Float64)

    return (1-(r^2)^(7.5))^7
end