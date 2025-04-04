
using MID
using LinearAlgebra

met =MID.MetT(zeros(3, 3), zeros(3, 3), zeros(3, 3, 3), zeros(3, 3, 3), 0.0, zeros(3))
B = MID.BFieldT(zeros(3), zeros(3), zeros(3, 3), zeros(3, 3), 0.0, zeros(3))
B_old = MID.BFieldT(zeros(3), zeros(3), zeros(3, 3), zeros(3, 3), 0.0, zeros(3))

rt = 0.1
θt = 0.2
R0 = 10.0
ζt = 0.3

MID.toroidal_metric!(met, rt, θt, ζt, R0);

geo = GeoParamsT(R0=R0)

isl = IslandT(m0=1, n0=2, A=0.01)
prob = init_problem(q=Axel_q, geo=geo)#, isl=isl); 

MID.MagneticField.compute_B!(B, met, prob.q, prob.isl, rt, θt, ζt)
MID.MagneticField.old_compute_B!(B_old, met, prob.q, prob.isl, rt, θt, ζt)


display(B.dB)
display(B_old.dB)
Γ = zeros(3, 3)

Γ[1, 1] = 1
Γ[2, 2] = 1
Γ[3, 3] = 1
    #display(Γ)


for i in 1:3
    #this is better but a fkn garbage way to do this!
    Γ[1, 1] -= met.gl[1, i] * B.b[i] * B.b[1]
    Γ[1, 2] -= met.gl[1, i] * B.b[i] * B.b[2]
    Γ[1, 3] -= met.gl[1, i] * B.b[i] * B.b[3]

    Γ[2, 1] -= met.gl[2, i] * B.b[i] * B.b[1]
    Γ[2, 2] -= met.gl[2, i] * B.b[i] * B.b[2]
    Γ[2, 3] -= met.gl[2, i] * B.b[i] * B.b[3]

    Γ[3, 1] -= met.gl[3, i] * B.b[i] * B.b[1]
    Γ[3, 2] -= met.gl[3, i] * B.b[i] * B.b[2]
    Γ[3, 3] -= met.gl[3, i] * B.b[i] * B.b[3]
end
display(Γ)

Γ = zeros(3, 3);
dΓ = zeros(3, 3, 3); #last ind is deriv.

for i in 1:3, j in 1:3
    Γ[i, j] = dδ(i, j) - dot(met.gl[i, :], B.b[:])*B.b[j]
end

display(Γ)


function dδ(i, j)
    if i==j
        return 1.0
    else
        return 0.0
    end
end

#so the new W gives some complete gib.
#not sure what is going on tbh.
old = MID.WeakForm.compute_Tj(met, B)

old_combo = zeros(9, 9)

for μ=1:9, i=1:3
    old_combo[i, μ] += old[i, μ] 
    old_combo[μ, i] += old[i, μ] 
end 

display(old_combo)

#while different, seems to be converging to the same damping value as everything else
new = MID.WeakForm.new_compute_Tj(met, B)

#added some islands for extra spice, diff is all over the place now.
#still often very small, yet there are some cases which are exactly zero, seems inconsistent???
#perhaps the lc tensor has some v small aspect to it? No evidence of this though!
#diff is indep of metric???? seems odd? maybe it is just floating point error?
diff = new .- old
println(diff)

#maybe we need to combine our methods and expand the lc and the other contractions, ie just write big ol sums... will be a nightmare though!
#other option is to start comparing cases with Mathematica, probably easiest to consider cylindrical limit...


jparold = MID.WeakForm.jparonB(met, B)
jparnew = MID.WeakForm.new_jparonB(met, B)

jpar_diff = jparold-jparnew


Wold = zeros(ComplexF64, 9, 9, 1, 1, 1);

@views MID.WeakForm.compute_W!(Wold[:, :, 1, 1, 1], met, B)

Wnew = zeros(ComplexF64, 9, 9, 1, 1, 1);

@views MID.WeakForm.new_compute_W!(Wnew[:, :, 1, 1, 1], met, B)


display(Wnew .- Wold)




#so now Tl is identical across different cases.
old = MID.WeakForm.compute_Tl(met, B)

new = MID.WeakForm.new_compute_Tl(met, B)

diff = new .- old

diff[5, 6]
println(diff)


Iold = zeros(ComplexF64, 9, 9, 1, 1, 1);

Hold = @views MID.WeakForm.compute_I!(Iold[:, :, 1, 1, 1], met, B, 1.0, 1.0e-7)

Inew = zeros(ComplexF64, 9, 9, 1, 1, 1);

Hnew = @views MID.WeakForm.new_compute_I!(Inew[:, :, 1, 1, 1], met, B, 1.0, 1.0e-7)

#println(Hold - Hnew)

#this is extremly different! especially the complex value??
diff = Iold[:, :, 1, 1, 1] - Inew[:, :, 1, 1, 1]


out1 = zeros(9, 9)

for i in 1:9, j in 1:9
    out1[i, j] = Hnew[i] * Hnew[j]
end

display(out1)

out2 = Hnew * Hnew'

display(out2-out1)
