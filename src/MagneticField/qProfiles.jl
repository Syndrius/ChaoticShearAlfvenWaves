

"""

Not sure on the scope of this, as we want this to just be part of Magnetic field essentially, although it is useful for plotting as well.

#this file needs to be fixed up once verification has been done.
#currently a mess of different q's that we have tried to compare against.

"""

#removed islands from inputs now! Island data type has changed
function Axel_q(r::Float64) :: Tuple{Float64, Float64}

    q = @. 1.05 + 0.55 * r^2
    dq = @. 2*0.55 * r

    return q, dq
end


function comparison_bowden_q(r::Float64)
    q = 1.5+(2-1.5)*r^2
    dq = r

    return q, dq

end

function singular_bowden_q(r::Float64)

    q = @. 1+(3-1)*r^2
    dq = @. 4*r

    return q, dq
end

function fu_dam_q(r::Float64)
    q = @. 1 + r^2
    dq = @. 2*r
    return q, dq
end

#this is way better, passing isl into the q-profile seems to cook this massivly.
#may just have to come up with a few q-profiles.
#still gives runtime dispatch warnings, but doesn't take up all of the time.
function default_island_q(r::Float64)

    q = 1 / (1 / (5/2) - 4 / (5/2)^2 * (r^2/2-0.125))
    dq = 4*(5/2)^2* 4 * r / (2*(5/2)-4*r^2 + 2*4*0.125)^2

    return q, dq
end