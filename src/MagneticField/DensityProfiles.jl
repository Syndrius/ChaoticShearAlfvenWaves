

#not really sure where this should go, but it is closest to q-profiles!
#mainly just used for verifying damping rates.

function axel_dens(r)
    #Axel's density function
    return 1/2 * (1 - tanh((r-0.8)/0.1))

end

function uniform_dens(r)
    return 1.0
end

function bowden_singular_dens(r)
    return 1/2 * (1-tanh((r-0.7)/0.05))
end

function comparison_bowden_dens(r)

    return (1-(r^2)^(7.5))^7
end