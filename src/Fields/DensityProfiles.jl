"""
Profile allows a GAE to form when combined with gae_q.
Taken from Van Rij et al. 1985.
"""
function gae_dens(r::Float64)

    return 1-r^2/(1+0.1)
end


"""
    uniform_dens(x1::Float64)

Uniformd density function, default density, returns 1.0.
"""
function uniform_dens(x1::Float64)
    return 1.0
end


"""
Density profile for comparison for continuum damping, when paired with damping_q.
Profile taken from Bowden and Hole 2015.
"""
function damping_dens(x1::Float64)
    return 1/2 * (1-tanh((x1-0.7)/0.05))
end

