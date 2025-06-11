#testing that the construction of I is doing exactly what we expect. Comparison with Mathematica.

using MID
using MID.WeakForm
using MID.Geometry
using MID.Equilibrium
using LinearAlgebra
#%%

function test_flr_I(r, θ, ζ, prob)

    I = zeros(ComplexF64, 9, 9)
    met = MetT()
    B = BFieldT()
    #tm = TM()
    F = zeros(9)
    toroidal_metric!(met, r, θ, ζ, prob.geo.R0)

    compute_B!(B, met, prob.q, prob.isls, r, θ, ζ)

   
    MID.WeakForm.compute_F!(B, met, F)
    #views are cooking julia
    #MID.WeakForm.compute_I!(@view I[:, :, 1, 1, 1],  met, B, 1.0, prob.flr, tm.D, tm.F)

    δ = 0.01

    n = 1.0
    for j=1:9, i=1:9
        #computes the non-ideal effects
        #assumes that Te ≈ Ti
        #then ρ_s = ρ_i
        #δ is artifical resitivity
        #ρ_i is ion gyro radius
        #δ_e is electron resistivity
        I[i, j] = -F[i] * F[j] * met.J[1] * n / B.mag_B[1]^2*(1.0im * δ )
    end

    return I
end
#this seems to be fine.
function test_I(r, θ, ζ, prob)

    I = zeros(ComplexF64, 3, 3)
    met = MetT()
    B = BFieldT()
    tm = TM()
    toroidal_metric!(met, r, θ, ζ, prob.geo.R0)

    compute_B!(B, met, prob.q, prob.isls, r, θ, ζ)
   
    MID.WeakForm.compute_D!(B, met, tm.D)
    #views are cooking julia
    #MID.WeakForm.compute_I!(@view I[:, :, 1, 1, 1],  met, B, 1.0, prob.flr, tm.D, tm.F)

    n = 1.0
    I[1:3, 1:3] = tm.D .* (n * met.J[1] / B.mag_B[1]^2)

    return I

end

#looks ok, assuming the weakform is actually correct.
function test_Tl(r, θ, ζ, prob)

    I = zeros(ComplexF64, 3, 3, 1, 1, 1)
    W = zeros(ComplexF64, 9, 9)
    met = MetT()
    B = BFieldT()
    tm = TM()
    toroidal_metric!(met, r, θ, ζ, prob.geo.R0)

    compute_B!(B, met, prob.q, prob.isls, r, θ, ζ)
   
    MID.WeakForm.compute_D!(B, met, tm.D)
    #@views MID.WeakForm.compute_I!(I[:, :, 1, 1, 1],  met, B, 1.0, prob.flr, tm.D, tm.F)

    MID.WeakForm.compute_C!(B, tm.C)
    #compute_C_old!(B, C)

    mul!(tm.T, tm.C', tm.D)

    mul!(W, tm.T, tm.C)

    W .*= met.J * B.mag_B[1]^2

    return W

end
#%%
#so this is fine, at least numerically. As long as we trust the thesis. this expression actually seems fine.
#wild tbh.
function test_Tj(r, θ, ζ, prob)

    W = zeros(ComplexF64, 9, 9)
    met = MetT()
    B = BFieldT()
    tm = TM()
    toroidal_metric!(met, r, θ, ζ, prob.geo.R0)

    compute_B!(B, met, prob.q, prob.isls, r, θ, ζ)
   
    sf = - met.J[1] * MID.WeakForm.jparonB(B, met) / 2.0

    Γ = tm.Γ
    dΓ = tm.dΓ
    K = tm.K

    lct = MID.WeakForm.get_lc_tensor()
    #computes the Γ matrix and its derivative dΓ
    MID.WeakForm.compute_Γ!(B, met, Γ, dΓ)

    #display("Gam is")
    #display(Γ)
    val = 0

    for n in 1:3

        #make sure K is zero each time.
        #alternative would be to set K values in first iteration of
        #loop, then have a smaller loop.
        fill!(K, 0.0)
        #computes K, contribution to W for ∂^2 terms without Γ derivatives.
        MID.WeakForm.compute_K!(met, Γ, K, n)

       

        #this loop computes the contributions for first derivatives,
        #so contains the derivatives of Γ.
        for q in 1:3
            #ideally this will have a better and more meaningful name.
            val = 0
            for i in 1:3, j in 1:3, k in 1:3
                val += Γ[i, n] / met.J[1] * lct[i, j, k] * dΓ[k, q, j]
            end

            #transpose is added to reflect that the two terms of Tj are identical expect Ψ -> Φ.
            W[n, q] += val * sf
            W[q, n] += val * sf
        end

        
        #add the K contributions.
        #transpose is added to reflect that the two terms of Tj are identical expect Ψ -> Φ.
        W[n, 4:9] .+= K .* sf
        W[4:9, n] .+= K .* sf

        #display("Wloc part is")
        #println(W[4:6, n])
        #println(W[7:9, n])

    end
    return W
end

#%%



geo = init_geo(R0=10.0)
isl = init_island(m0=3, n0=-2, A=0.0001)
prob = init_problem(q = fu_dam_q, isl=isl, geo=geo)
#%%
rtest = 0.34
θtest = 1.45
ζtest = 4.53

#so assuming (big assume!) our maths is actually correct and our equations are valid
#by comparison with mathematica, we can be pretty confident that the numerical implementation has been done properly.
res = test_I(rtest, θtest, ζtest, prob)
res = test_Tl(rtest, θtest, ζtest, prob)
res = test_Tj(rtest, θtest, ζtest, prob)
res = test_flr_I(rtest, θtest, ζtest, prob)
display(imag.(res))
