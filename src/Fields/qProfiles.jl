"""
Simple quadratic q-profile.
"""
function quadratic_q(x1::Float64)
    q = @. 1 + x1^2
    dq = @. 2*x1
    return q, dq
end

"""
Designed to mimic a normal q-profile, q(0)=1, q(1)=2
while making the iota profile linear in ψ
this has 3/2 island at ψ=0.5, and 4/3 island at ψ=2/3
This q-profile is used for QFM and periodic orbits, as the algorithm is simplified.
"""
function cantori_q(ψ::Float64)
    return 1/(1-0.5*ψ), 1/2 / (1-1/2*ψ)^2
end


"""
Profile to make analytical coordinates for inside an magnetic island.
"""
function island_q(κ::Float64, isl::CoordIslandT)
    K, E = Elliptic.ellipke(κ)

    q = isl.w / (isl.m0*π) * K

    dq = isl.w / (isl.m0*π) * (E - (1-κ) * K) / (2*(1-κ)*κ)

    return q, dq
end

"""
Profile to match the analytical q-profile for islands, for geometric radius.
"""
function island_q(r::Float64, isl::RadialIslandT)
    q0 = isl.q0
    qp = isl.qp 
    r0 = isl.r0
    q = 1 / (1 / q0 - qp / (2*q0^2*r0) * (r^2-r0^2))
    dq = 4*qp*q0^2*r*r0 / (2*q0*r0 - qp * (r^2-r0^2))^2
    return q, dq
end

"""
Profile to match the analytical q-profile for islands, for toroidal flux.
"""
function island_q(ψ::Float64, isl::FluxIslandT)
    q0 = isl.q0
    qp = isl.qp 
    ψ0 = isl.ψ0
    q = 1 / (1/q0 - qp/q0^2*(ψ-ψ0))
    dq = qp/q0^2 * q^2
    return q, dq
end


"""
Profile allows a GAE to form when combined with gae_dens.
Taken from Van Rij et al. 1985.
"""
function gae_q(r::Float64)
    q = 1 + 2*r^2
    dq = 4*r
    return q, dq
end


"""
Q profile for comparison for continuum damping, when paired with damping_dens.
Profile taken from Bowden and Hole 2015.
"""
function damping_q(x1::Float64)

    q = 1+(3-1)*x1^2
    dq = 4*x1

    return q, dq
end

