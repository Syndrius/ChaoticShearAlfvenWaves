

"""
Devised q-profile such that gap and 5/4 island both occur at r=0.5.

# Args
r::Float64 - Radial point.
"""
function island_damping_q(r)
    a = 1.15
    b = 0.4
    q = a + b * r^2
    dq = 2 * b * r
    return q, dq
end

"""
Q profile from Axel's paper, should be renamed.

# Args
r::Float64 - Radial point.
"""
function Axel_q(r::Float64) 

    q = @. 1.05 + 0.55 * r^2
    dq = @. 2*0.55 * r

    return q, dq
end

"""
Q profile from Bowden comparison paper, should be renamed, or removed.

# Args
r::Float64 - Radial point.
"""
function comparison_bowden_q(r::Float64)
    q = 1.5+(2-1.5)*r^2
    dq = r

    return q, dq

end

"""
Q profile from Bowden singular paper, should be renamed.

# Args
r::Float64 - Radial point.
"""
function singular_bowden_q(r::Float64)

    q = 1+(3-1)*r^2
    dq = 4*r

    return q, dq
end

"""
Q profile from Fu and Van Dam paper, should be renamed.

# Args
r::Float64 - Radial point.
"""
function fu_dam_q(r::Float64)
    q = @. 1 + r^2
    dq = @. 2*r
    return q, dq
end

"""
Q profile from Qu's paper, should be renamed.

# Args
r::Float64 - Radial point.
"""
function default_island_q(r::Float64)

    q = 1 / (1 / (5/2) - 4 / (5/2)^2 * (r^2/2-0.125))
    dq = 4*(5/2)^2* 4 * r / (2*(5/2)-4*r^2 + 2*4*0.125)^2

    return q, dq
end