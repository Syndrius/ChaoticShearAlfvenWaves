
#testing the isalnd amplitude function to prevent siliness at r=0
#main idea is to make sure B^r << B^θ, B^ζ

#i.e. given some q-profile, in this case island_mode_21, what is the behaviour of each component (physical and conravarient) as r-> 0

N = 100

ramp = LinRange(0, 1, N)

q = zeros(N)
for (i, r) in enumerate(ramp)
    #display(i)
    q[i], _ = island_mode_21(r)
end

plot(ramp, q)

#start with covarient parts
#assuming cylindrical geom...

#this seems to be an ok choice?? I think???

function isl_amp(r)
    if r > 0.3 #this will assume the island is at r=0.5
        return 1
    elseif r > 0.2
        #function chosen to be zero at r=0.2, 1 at r=0.3, and gradient of 0 at r=0.3.
        #subject to change obvs.
        return -100*r^2 + 60*r - 8
    else

        return 0
    end
end

Br = zeros(N)
Bθ = zeros(N)
Bζ = zeros(N)
R0 = 10.0
A = 1e-4
for (i, r) in enumerate(ramp)
    J = r * R0 #not sure what R0 should be...
    Br[i] = 1 / J * A * isl_amp(r) * 2 * 1
    Bθ[i] = r / (J * q[i])
    Bζ[i] = r / J
end

#ok so it is probably not Br vs Bθ, probably more A->0 faster than J -> zero.
#so given J~r, need A~r^2 I think.

plot(ramp, Br)
plot(ramp, Bθ)
plot(ramp, Bζ)


isl = IslandT(m0=2, n0=-1, w=0.03)

prob = init_problem(q=island_mode_21, met=cylindrical_metric!, geo=GeoParamsT(R0=1000), isl=isl)

display(prob.isl.A)