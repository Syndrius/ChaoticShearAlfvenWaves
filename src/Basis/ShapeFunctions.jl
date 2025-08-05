"""
Cubic Hermite splines consist of 4 shape functions, typically defined for t∈[0, 1].
To match the local domain of ξ∈[-1, 1] for Gaussian quadrature, we transform the shape functions to this domain.

The shape functions are defined as h_{ij}
where i=0, 1 dictates the basis, with 0 being the normal basis, and 1 being the tangent basis,
and j=0, 1 dicates the grid, meaning the important value is defined on the left edge (j=0) or right edge (j=1)

By important value, we not that each shape function has either value 1 or derivative 1 on one of the edges.

h00 has value 1 at ξ=-1, and is 0 at ξ=1, with derivative 0 at both edges.
h10 has derivative 1 at ξ=-1, zero for all other cases
h01 has value 1 at ξ=1
h11 has derivative 1 at ξ=1

Note that the tangent functions, h10 and h11, additionally scaling is required when transforming the domain,
this ensures the above behaviour occurs.
"""


#h00(-1) = 1
function h00(ξ)
    return (2 - 3*ξ + ξ^3) / 4
end

#h10'(-1) = 1
function h10(ξ)
    return (ξ - 1)^2 * (ξ + 1) / 4
end

#h01(1) = 1
function h01(ξ)
    return (2 + 3*ξ - ξ^3) / 4
end

#h11'(1) = 1
function h11(ξ)
    return (ξ - 1) * (ξ + 1)^2 / 4
end


#need to be more confident about this!
function dh00(ξ)

    return 3*(ξ^2 - 1) / 4
end

#h10'(-1) = 1
function dh10(ξ)
    return (3*ξ^2 - 2*ξ - 1) / 4
end

#h01(1) = 1
function dh01(ξ)
    return -3 * (ξ^2 - 1) / 4
end

#h11'(1) = 1
function dh11(ξ)
    return (3*ξ^2 + 2*ξ - 1) /4
end


function ddh00(ξ)

    return 3*ξ / 2
end

#h10'(-1) = 1
function ddh10(ξ)
    return (3*ξ - 1) / 2
end

#h01(1) = 1
function ddh01(ξ)
    return -3*ξ / 2
end

#h11'(1) = 1
function ddh11(ξ)
    return (1+3*ξ) / 2
end
