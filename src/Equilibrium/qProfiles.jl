#TODO -> should change these names and probably provide the function in the doc.@

#chosen so 1/1 surface is at r~0.05 and 5/2 surface is at r~0.95, should allow us to view the entire domain
#this results in a 4/3 island at r=0.45, and a 3/2 island at r=0.55.
function qfm_q(r::Float64)
    #a = 0.995833333
    #b = 1.66666666
    a = 239 / 240
    b = 5 / 3
    return a + b * r^2, 2 * b * r
end


#chosen to be 1+2*r^2, but shifted a bit so that the qfm surfaces at r=0 (1, 1) and r=1 (3, 1) do not go outside the domain.
#this is causing heaps of issues, unsure why. Seems like most of the problemo is coming from 
#outside the chaotic region, which is very weird.
#perhaps the islands in this case are too far away, meaning chaos only occurs when we have very large islands, causing problems in the outer regions.
#think perhaps our ideal q profile would be the (4, 3) and (3, 2) islands, starting at 1 ish
#but we probably need the two island chains to be within 0.1 of each other.
#this version is better, but we still have huge issues for r < 0.4
function old_qfm_q(r::Float64)
    a = 0.954545
    b = 1.515151
    return a + b * r^2, 2 * b * r
end

function chaos_q(r::Float64)

    return 1.0 + 2.0 * r^2, 4.0 * r
end

function qfm_benchmark_q(r::Float64)
    #a = 0.954545
    #b = 1.515151
    
    a = 1.93333
    b = 1.66666

    return a + b*r^2, 2 * b * r
end


function tae_isl_damping_q(r::Float64)
    #should give tae (4/5, -2) tae at 0.75 @~0.2 with (2, -1) isl at 0.5.
    a = 1.8
    b = 0.8
    q = a + b * r^2
    dq = 2 * b * r
    return q, dq
end

#IslandT(2, -1, 5.625e-5, 2.0, 2.0, 0.5, 0.03)
function inside_island_q(r)
    #A = 0.00015625000000000003
    #w = 0.05
    w = 0.03
    A = 5.625e-5
    m0 = 2
    return -w/(2*A*π*m0) * Elliptic.K(r), 0
end

#island mode case for 21 island.
function island_mode_21(r::Float64)
    q0 = 2/1
    qp = 2 #chosen pretty arbitrarily based on vibes of continuum.
    r0 = 0.5
    #ψ0 = 0.125
    #q = 1 / (1 / q0 - qp / (q0)^2 * (r^2/2-ψ0))
    #dq = 4*(q0)^2* qp * r / (2*(q0)-qp*r^2 + 2*qp*ψ0)^2
    q = 1 / (1 / q0 - qp / (2*q0^2*r0) * (r^2-r0^2))
    dq = 4*qp*q0^2*r*r0 / (2*q0*r0 - qp * (r^2-r0^2))^2
    return q, dq
end


#very similar to island_mode_21 but with modified qp for the gae.
function gae_isl_q(r::Float64)
    q0 = 2/1
    qp = 2 #chosen pretty arbitrarily based on vibes of continuum.
    qp = -8
    r0 = 0.5
    #ψ0 = 0.125
    #q = 1 / (1 / q0 - qp / (q0)^2 * (r^2/2-ψ0))
    #dq = 4*(q0)^2* qp * r / (2*(q0)-qp*r^2 + 2*qp*ψ0)^2
    q = 1 / (1 / q0 - qp / (2*q0^2*r0) * (r^2-r0^2))
    dq = 4*qp*q0^2*r*r0 / (2*q0*r0 - qp * (r^2-r0^2))^2
    return q, dq
end


#q profile based of Axel's 2024 iota profile
#similar to Zhisongs case, but we now have a 5/5 island.
#island is located at ~0.859
#ideally the island would be a wee bit closer to r=0.5...
function Axel_island_q(r::Float64)
    #note 0.7381 is r0^2/a^2, 
    #our case has a=1, not true for Axel...
    r02 = 0.7381
    ιp = 0.2039716
    q = 1 / (1 + ιp * (r^2 - r02))
    dq = -2 * ιp * r / (ιp * (r^2-r02) + 1)^2

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


function generic_island_q(r::Float64, isl::IslandT)

    q0 = isl.q0
    qp = isl.qp
    r0 = isl.r0
    #need to verify this!
    #maybe this will fix the flux surface problemo.
    q = 1 / (1 / q0 - qp / (2*q0^2*r0) * (r^2-r0^2))
    dq = 4*qp*q0^2*r*r0 / (2*q0*r0 - qp * (r^2-r0^2))^2
    return q, dq


end

#based on Zhisongs q-profile, but modifed to have m0=3 rather than 5 
#to increase the resolution of island modes.
function island_mode_q(r::Float64)
    q0 = 3/2
    qp = 1.6
    ψ0 = 0.125
    r0 = sqrt(2*ψ0)
    #q = 1 / (1 / q0 - qp / (q0)^2 * (r^2/2-ψ0))
    #dq = 4*(q0)^2* qp * r / (2*(q0)-qp*r^2 + 2*qp*ψ0)^2
    q = 1 / (1 / q0 - qp / (2*q0^2*r0) * (r^2-r0^2))
    dq = 4*qp*q0^2*r*r0 / (2*q0*r0 - qp * (r^2-r0^2))^2
    return q, dq
end



"""
    island_damping_q(r)

Devised q-profile such that gap and 5/4 island both occur at r=0.5.
"""
function island_damping_q(r)
    a = 1.15
    b = 0.4
    q = a + b * r^2
    dq = 2 * b * r
    return q, dq
end


function flr_q(r)
    a = 1.05
    b = 0.55

    q = a + b*r^2
    dq = 2 * b * r
    return q, dq

end


function symmetric_q(r)
    a = 1
    q = 1 + a * r * (1-r)
    dq = a - 2*a * r
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

    q = a+b*r^2
    dq = 2*b*r
    return q, dq
end


function test_q(r)

    #with 6/5 and 8/15 (3, 2) island will form at r=0.75,
    #and (5/6, 4) tae will form at r~0.57
    a = 6/5
    b = 8/15
    q = a+b*r^2
    dq = 2 * b * r
    return q, dq
end


"""
    Axel_q(r::Float64) 

Q profile from Axel's paper, should be renamed.
"""
function Axel_q(r::Float64) 

    q = @. 1.05 + 0.55 * r^2
    dq = @. 2*0.55 * r

    return q, dq
end


"""
    Axel_q(r::Float64) 

Q profile from Axel's paper, should be renamed.
"""
function flux_Axel_q(ψ::Float64) 

    q = @. 1.05 + 0.55 * 2 * ψ
    dq = @. 2*0.55 

    return q, dq
end

"""
    comparison_bowden_q(r::Float64)

Q profile from Bowden comparison paper, should be renamed, or removed.
"""
function comparison_bowden_q(r::Float64)
    q = 1.5+(2-1.5)*r^2
    dq = r

    return q, dq

end

"""
    bowden_singular_q(r::Float64)

Q profile from Bowden singular paper, should be renamed.
"""
function bowden_singular_q(r::Float64)

    q = 1+(3-1)*r^2
    dq = 4*r

    return q, dq
end

"""
    fu_dam_q(r::Float64)

Q profile from Fu and Van Dam paper, should be renamed.
"""
function fu_dam_q(r::Float64)
    q = @. 1 + r^2
    dq = @. 2*r
    return q, dq
end

"""
    default_island_q(r::Float64)

Q profile from Qu's paper, should be renamed.
"""
function default_island_q(r::Float64)

    q = 1 / (1 / (5/2) - 4 / (5/2)^2 * (r^2/2-0.125))
    dq = 4*(5/2)^2* 4 * r / (2*(5/2)-4*r^2 + 2*4*0.125)^2

    return q, dq
end
