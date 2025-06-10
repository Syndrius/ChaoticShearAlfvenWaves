#TODO -> should change these names and probably provide the function in the doc.@

#q profile required to use island coordinates
function island_coords_q(κ::Float64, isl::IslandT)

    K, E = Elliptic.ellipke(κ)

    q = isl.w / (isl.m0*π) * K

    dq = isl.w / (isl.m0*π) * (E - (1-κ) * K) / (2*(1-κ)*κ)

    return q, dq
end

#this is a terrible solution but oh wel
function island_coords_21a_q(κ::Float64)

    w = 0.03
    m0 = 2
    
    K, E = Elliptic.ellipke(κ)

    q = w / (m0*π) * K

    dq = w / (m0*π) * (E - (1-κ) * K) / (2*(1-κ)*κ)

    return q, dq
end


#Toroidal equivalent ot eh island coords q, allows for direct comparison.
function island_equiv_q(r::Float64, isl::IslandT)

    q0 = isl.q0
    qp = isl.qp #chosen pretty arbitrarily based on vibes of continuum.
    r0 = isl.r0
    #this form is to allows analytical integration in the construction of island coordinates.
    q = 1 / (1 / q0 - qp / (2*q0^2*r0) * (r^2-r0^2))
    dq = 4*qp*q0^2*r*r0 / (2*q0*r0 - qp * (r^2-r0^2))^2
    return q, dq
end

#JLD2 cannot save anonymouse functions, so we probably need to hard code them unfort
function island_21_q(r::Float64)
    q0 = 2/1
    qp = 2.0
    r0 = 0.5
    #this form is to allows analytical integration in the construction of island coordinates.
    q = 1 / (1 / q0 - qp / (2*q0^2*r0) * (r^2-r0^2))
    dq = 4*qp*q0^2*r*r0 / (2*q0*r0 - qp * (r^2-r0^2))^2
    return q, dq
end

function island_32_q(r::Float64)
    q0 = 3/2
    qp = 2.0
    r0 = 0.5
    #this form is to allows analytical integration in the construction of island coordinates.
    q = 1 / (1 / q0 - qp / (2*q0^2*r0) * (r^2-r0^2))
    dq = 4*qp*q0^2*r*r0 / (2*q0*r0 - qp * (r^2-r0^2))^2
    return q, dq
end
#
#desgined for chaotic region between ~0.77-0.9
#with (7, -4), (5, -3), (8, -5) islands
function low_shear_qfm_q(x1::Float64)
    a = 1.2
    b = 0.6
    c = 0.1
    return a + b * x1^2 + c*x1^4, 2 * b * x1 + 4*c*x1^3
end

#stupid af q-profile, should have n=3 islands, starting from m=7 at 0.7
#up to m=14 at 0.9
function single_n_q_prof(x1::Float64)
    a = 1.05
    b = 1.3
    c = 0.83254
    d = 3.8
    return a + b*x1^2 + c*x1^4 + d*x1^6, 2 * b * x1 + 4*c*x1^3 + 6*d*x1^5
    #return a + b*x1^2 + c*x1^10 + d*x1^4, 2*b*r + 10*c*x1^9 + 4*d*x1^3
    #return a + b*x1^6, 6 * b * x1^5
    #return a + b*x1^4, 4 * b * x1^3
end

function low_shear_q(x1::Float64)
    a = 1.05
    b = 0.45
    return a + b * x1^2, 2 * b * x1
end

#chosen so that 4, 3 island at r=0.4, (3, 2) island at 0.6, (1, 1)/(2, 1) tae at r~0.5
#this gives a=6/5 and b=5/6
#q profile is then shifted so that end points are moved to 0.05 and 0.95 for qfm surfaces
function qfm_q(x1::Float64)
    a = 539/450
    b = 8/9
    return a + b * x1^2, 2 * b * x1
end

#chosen so 1/1 surface is at r~0.05 and 5/2 surface is at r~0.95, should allow us to view the entire domain
#this results in a 4/3 island at r=0.45, and a 3/2 island at r=0.55.
function old_qfm_q(x1::Float64)
    #a = 0.995833333
    #b = 1.66666666
    a = 239 / 240
    b = 5 / 3
    return a + b * x1^2, 2 * b * x1
end


#chosen to be 1+2*x1^2, but shifted a bit so that the qfm surfaces at r=0 (1, 1) and r=1 (3, 1) do not go outside the domain.
#this is causing heaps of issues, unsure why. Seems like most of the problemo is coming from 
#outside the chaotic region, which is very weird.
#perhaps the islands in this case are too far away, meaning chaos only occurs when we have very large islands, causing problems in the outer regions.
#think perhaps our ideal q profile would be the (4, 3) and (3, 2) islands, starting at 1 ish
#but we probably need the two island chains to be within 0.1 of each other.
#this version is better, but we still have huge issues for r < 0.4
function oldest_qfm_q(x1::Float64)
    a = 0.954545
    b = 1.515151
    return a + b * x1^2, 2 * b * x1
end

function chaos_q(x1::Float64)

    return 1.0 + 2.0 * x1^2, 4.0 * x1
end

function qfm_benchmark_q(x1::Float64)
    #a = 0.954545
    #b = 1.515151
    
    a = 1.93333
    b = 1.66666

    return a + b*x1^2, 2 * b * x1
end


function tae_isl_damping_q(x1::Float64)
    #should give tae (4/5, -2) tae at 0.75 @~0.2 with (2, -1) isl at 0.5.
    a = 1.8
    b = 0.8
    q = a + b * x1^2
    dq = 2 * b * x1
    return q, dq
end

#IslandT(2, -1, 5.625e-5, 2.0, 2.0, 0.5, 0.03)
function inside_island_q(r)
    #A = 0.00015625000000000003
    #w = 0.05
    w = 0.03
    A = 5.625e-5
    m0 = 2
    return -w/(2*A*π*m0) * Elliptic.K(x1), 0
end

#island mode case for 21 island.
function island_mode_21(x1::Float64)
    q0 = 2/1
    qp = 2 #chosen pretty arbitrarily based on vibes of continuum.
    r0 = 0.5 #this may still need to change to x10, but we are actually using r for this case
    #ψ0 = 0.125
    #q = 1 / (1 / q0 - qp / (q0)^2 * (x1^2/2-ψ0))
    #dq = 4*(q0)^2* qp * x1 / (2*(q0)-qp*x1^2 + 2*qp*ψ0)^2
    q = 1 / (1 / q0 - qp / (2*q0^2*r0) * (x1^2-r0^2))
    dq = 4*qp*q0^2*x1*r0 / (2*q0*r0 - qp * (x1^2-r0^2))^2
    return q, dq
end


#very similar to island_mode_21 but with modified qp for the gae.
function gae_isl_q(x1::Float64)
    q0 = 2/1
    qp = 2 #chosen pretty arbitrarily based on vibes of continuum.
    qp = -8
    r0 = 0.5
    #ψ0 = 0.125
    #q = 1 / (1 / q0 - qp / (q0)^2 * (x1^2/2-ψ0))
    #dq = 4*(q0)^2* qp * x1 / (2*(q0)-qp*x1^2 + 2*qp*ψ0)^2
    q = 1 / (1 / q0 - qp / (2*q0^2*r0) * (x1^2-r0^2))
    dq = 4*qp*q0^2*r*r0 / (2*q0*r0 - qp * (x1^2-r0^2))^2
    return q, dq
end


#q profile based of Axel's 2024 iota profile
#similar to Zhisongs case, but we now have a 5/5 island.
#island is located at ~0.859
#ideally the island would be a wee bit closer to r=0.5...
function Axel_island_q(x1::Float64)
    #note 0.7381 is r0^2/a^2, 
    #our case has a=1, not true for Axel...
    r02 = 0.7381
    ιp = 0.2039716
    q = 1 / (1 + ιp * (x1^2 - r02))
    dq = -2 * ιp * x1 / (ιp * (x1^2-r02) + 1)^2

    return q, dq
end

function flux_fu_dam_q(ψ::Float64)
    q = @. 1 + ψ * 2
    dq = @. 2
    return q, dq
end


function flux_generic_island_q(ψ::Float64, isl::IslandT)

    q0 = isl.q0
    qp = isl.qp
    ψ0 = isl.ψ0
    #need to verify this!
    q = 1 / (1 / q0 - qp / (q0)^2 * (ψ-ψ0))
    dq = (q0)^2 * qp / (q0 - qp*(ψ - ψ0))^2
    return q, dq


end


function generic_island_q(x1::Float64, isl::IslandT)

    q0 = isl.q0
    qp = isl.qp
    r0 = isl.r0
    #need to verify this!
    #maybe this will fix the flux surface problemo.
    q = 1 / (1 / q0 - qp / (2*q0^2*r0) * (x1^2-r0^2))
    dq = 4*qp*q0^2*r*r0 / (2*q0*r0 - qp * (x1^2-r0^2))^2
    return q, dq


end

#based on Zhisongs q-profile, but modifed to have m0=3 rather than 5 
#to increase the resolution of island modes.
function island_mode_q(x1::Float64)
    q0 = 3/2
    qp = 1.6
    ψ0 = 0.125
    r0 = sqrt(2*ψ0)
    #q = 1 / (1 / q0 - qp / (q0)^2 * (x1^2/2-ψ0))
    #dq = 4*(q0)^2* qp * x1 / (2*(q0)-qp*x1^2 + 2*qp*ψ0)^2
    q = 1 / (1 / q0 - qp / (2*q0^2*r0) * (x1^2-r0^2))
    dq = 4*qp*q0^2*r*r0 / (2*q0*r0 - qp * (x1^2-r0^2))^2
    return q, dq
end



"""
    island_damping_q(r)

Devised q-profile such that gap and 5/4 island both occur at r=0.5.
"""
function island_damping_q(r)
    a = 1.15
    b = 0.4
    q = a + b * x1^2
    dq = 2 * b * x1
    return q, dq
end


function flr_q(r)
    a = 1.05
    b = 0.55

    q = a + b*x1^2
    dq = 2 * b * x1
    return q, dq

end


function symmetric_q(r)
    a = 1
    q = 1 + a * x1 * (1-r)
    dq = a - 2*a * x1
    return q, dq
end


#q-profile designed for tae at r=0.5 for (2, -2), (3, -2) with a (3,2) island at r=0.75.
#used to be 
#a=1.05
#b=0.8
#changed to form such that m0=n of tae, to try and get overlap.
function island_3_2_q(r)
    a = 1.4
    b = 0.4

    q = a+b*x1^2
    dq = 2*b*r
    return q, dq
end


function test_q(r)

    #with 6/5 and 8/15 (3, 2) island will form at r=0.75,
    #and (5/6, 4) tae will form at r~0.57
    a = 6/5
    b = 8/15
    q = a+b*x1^2
    dq = 2 * b * x1
    return q, dq
end


"""
    Axel_q(x1::Float64) 

Q profile from Axel's paper, should be renamed.
"""
function Axel_q(x1::Float64) 

    q = @. 1.05 + 0.55 * x1^2
    dq = @. 2*0.55 * x1

    return q, dq
end


"""
    Axel_q(x1::Float64) 

Q profile from Axel's paper, should be renamed.
"""
function flux_Axel_q(ψ::Float64) 

    q = @. 1.05 + 0.55 * 2 * ψ
    dq = @. 2*0.55 

    return q, dq
end

"""
    comparison_bowden_q(x1::Float64)

Q profile from Bowden comparison paper, should be renamed, or removed.
"""
function comparison_bowden_q(x1::Float64)
    q = 1.5+(2-1.5)*x1^2
    dq = r

    return q, dq

end

"""
    bowden_singular_q(x1::Float64)

Q profile from Bowden singular paper, should be renamed.
"""
function bowden_singular_q(x1::Float64)

    q = 1+(3-1)*x1^2
    dq = 4*r

    return q, dq
end

"""
    fu_dam_q(x1::Float64)

Q profile from Fu and Van Dam paper, should be renamed.
"""
function fu_dam_q(x1::Float64)
    q = @. 1 + x1^2
    dq = @. 2*x1
    return q, dq
end

"""
    default_island_q(x1::Float64)

Q profile from Qu's paper, should be renamed.
"""
function default_island_q(x1::Float64)

    q = 1 / (1 / (5/2) - 4 / (5/2)^2 * (x1^2/2-0.125))
    dq = 4*(5/2)^2* 4 * x1 / (2*(5/2)-4*x1^2 + 2*4*0.125)^2

    return q, dq
end
