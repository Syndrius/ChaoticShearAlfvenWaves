

#seeif we can replicate some of the literature results using Berks equation
#and see if we can determine that the difference in equation is really the problem.

#code is taken from fu-dam stuff earlier, need to add damping and change this to match our base case
#think Bowden singular is the best test example as it gives exact tae freq.
#but maybe axel will be better as he avoids the m=0 case.
#start with Axel_q for easiest boundary conditions!
#we will also need to add in the desnity profile.

#probably want to split this up, and change to sparse for speed, takes ages atm.
#hoepfully we can replicate George's result.

using Plots; plotlyjs()
using FastGaussQuadrature
using MID
using LinearAlgebra
using Arpack


function grid_to_ind_coupled(i, N, h, eqn)

    grid_id = [0, 0, 1, 1]

    basis_id = [0, 1, 0, 1]

    return 2 * (eqn-1)*N + 1 + 2*(i-1) + basis_id[h] + 2*grid_id[h]

end


function phi_fudam(N, R, m=1, n=-1)
    #this works as expected, tae freq is a tiny bit different, that is probably owing to the different approximations used, ie we have removed the extra term, but is very close.
    #this  version is also slightly closer to Berk's quoted 0.31.

    m1 = m+1

    rgrid = LinRange(0, 1, N)
    gp = 5
    ξ, wg = gausslegendre(gp)
    H, dH, ddH = MID.Misc.hermite_basis(ξ)

    W = zeros(4*N, 4*N)
    I = zeros(4*N, 4*N)

    ϕm = zeros(3, 4, gp)
    ϕm1 = zeros(3, 4, gp)

    ϕm[1, :, :] = H
    ϕm[2, :, :] = dH
    ϕm[3, :, :] = ddH

    #probably not needed.
    ϕm1[1, :, :] = H
    ϕm1[2, :, :] = dH
    ϕm1[3, :, :] = ddH

    for i in 1:N-1

        r, dr = MID.Misc.local_to_global(i, ξ, collect(rgrid))

        #q = @. 1.05 + 0.55 * r^2 #Axels q-profile, use m=2, n=-2

        q = @. 1.0 + 2.0*r^2 #Bowden Singular q, use m=1, n=-1

        km = @. 1/R*(m/q + n)

        km1 = @. 1/R*(m1/q + n)

        jac = dr/2

        ϵ = @. 5*r/(2*R)

        #cylindrical limit!
        #ϵ = zeros(length(r))


        #equation 1
        for trial in 1:4

            trial_ind_m = grid_to_ind_coupled(i, N, trial, 1)
            trial_ind_m1 = grid_to_ind_coupled(i, N, trial, 2)
            for test in 1:4
                test_ind_m = grid_to_ind_coupled(i, N, test, 1)
                test_ind_m1 = grid_to_ind_coupled(i, N, test , 2)
                for j in 1:gp

                    #going to assume v_A = 1?
                    #eq1 I
                    #I's are negative as I and W need to be opposite signs for eval problem
                    # I believe this is actually wrong, didn't account for km being a function of r...
                    #we may have to work with Em for now...
                    #Axels form is essentially this equation but with ϕ, it contains a term of dk^2/dr, which would be the extra term we have missed.
                    #ideally, sticking with ϕ would be better just for completness.
                    #this will need a lot more work to get working, will have to see how we go
                    #first step should be to use Em form, and solve Bowden's case to see if we can replicate damping!
                    I[test_ind_m, trial_ind_m] -= -jac * ϕm[2, test, j] * r[j] * ϕm[2, trial, j] * 1/jac^2 * wg[j]
                    I[test_ind_m, trial_ind_m] -= -jac * (m^2) * ϕm[1, test, j] / r[j] * ϕm[1, trial, j] * wg[j]
                    I[test_ind_m, trial_ind_m1] -= -jac * ϕm[2, test, j] * r[j] * ϵ[j] * ϕm1[2, trial, j] * 1/jac^2 * wg[j]

                    #optional extra term for I goes here.

                    #eq1 W
                    W[test_ind_m, trial_ind_m] += jac * ϕm[2, test, j] * r[j] * km[j]^2 * ϕm[2, trial, j] * 1/jac^2 * wg[j]
                    W[test_ind_m, trial_ind_m] += jac * (m^2) * ϕm[1, test, j] / r[j] * km[j]^2 * ϕm[1, trial, j] * wg[j]
                    


                    #eq2 I
                    I[test_ind_m1, trial_ind_m1] -= -jac * ϕm1[2, test, j] * r[j] * ϕm1[2, trial, j] * 1/jac^2 * wg[j]
                    I[test_ind_m1, trial_ind_m1] -= -jac * (m1^2) * ϕm1[1, test, j] / r[j] * ϕm1[1, trial, j] * wg[j]
                    I[test_ind_m1, trial_ind_m] -= -jac * ϕm1[2, test, j] * r[j] * ϵ[j] * ϕm[2, trial, j] * 1/jac^2 * wg[j]
                    
                    #option I term goes here

                    #eq2 W
                    W[test_ind_m1, trial_ind_m1] += jac * ϕm1[2, test, j] * r[j] * km1[j]^2 * ϕm1[2, trial, j] * 1/jac^2 * wg[j]
                    W[test_ind_m1, trial_ind_m1] += jac * (m1^2) * ϕm1[1, test, j] / r[j] * km1[j]^2 * ϕm1[1, trial, j] * wg[j]


                end
            end
        end
    end

    #apply the boundary conditions.
    if abs(m)==0
        bc1 = [2, 2*N-1]
    else
        bc1 = [1, 2*N-1]
    end
    if abs(m1)==0
        bc2 = [2*N+2, 4*N-1]
    else
        bc2 = [2*N+2, 4*N-1]
    end

    bcinds = vcat(bc1, bc2)

    for i in bcinds
        W[i, :] .= 0.0
        W[:, i] .= 0.0
        I[i, :] .= 0.0
        I[:, i] .= 0.0
    end

    for i in bcinds
        W[i, i] = 1.0
        I[i, i] = 1.0
    end
    ω2, evals = eigen(Hermitian(W), Hermitian(I))

    Emsol = evals[1:2:2*N, :]
    Em1sol = evals[2*N+1:2:end, :]
    return R .*sqrt.(ω2), Emsol, Em1sol


end

#grid_to_ind_coupled(2, 10, 1, 2)


#all this should need is a q-profile
function fudam(N, R, m, n, rgrid, δ)

    #Fu's paper is unclear, instead we will solve Berks form and
    #use that to benchamrk our more general code
    m1 = m+1

    

    #note we assume Berks form for the current term is correct
    #and we use the additional assumption of km = -km+1
    #to simplify the perturbation term.

    #rgrid = LinRange(0, 1, N)
    gp = 5
    ξ, wg = gausslegendre(gp)
    H, dH, ddH = MID.Misc.hermite_basis(ξ)

    W = zeros(ComplexF64, 4*N, 4*N)
    I = zeros(ComplexF64, 4*N, 4*N)

    Em = zeros(3, 4, gp)
    Em1 = zeros(3, 4, gp)

    Em[1, :, :] = H
    Em[2, :, :] = dH
    Em[3, :, :] = ddH

    #probably not needed.
    Em1[1, :, :] = H
    Em1[2, :, :] = dH
    Em1[3, :, :] = ddH

    

    for i in 1:N-1

        r, dr = MID.Misc.local_to_global(i, ξ, collect(rgrid))

        q = @. 1.0 + 2.0*r^2

        #qp = @. 2*r

        dens = @. 1/2*(1-tanh((r-0.7)/0.05))

        ddens = @. -10 * sech((r-0.07)/0.05)^2

        qpp = 2*ones(length(r))

        km = @. 1/R*(m/q + n)

        #kmp = @. -m*qp/(R*q^2)

        #kmpp = @. m * (2*qp^2 - q*qpp)/(R*q^3)

        km1 = @. 1/R*(m1/q + n)

        #km1p = @. -m1*qp/(R*q^2)

        #km1pp = @. m1 * (2*qp^2 - q*qpp)/(R*q^3)

        jac = dr/2

        ϵ = @. 5*r/(2*R)

        #cylindrical limit!
        #ϵ = zeros(length(r))


        #equation 1
        for trial in 1:4

            trial_ind_m = grid_to_ind_coupled(i, N, trial, 1)
            trial_ind_m1 = grid_to_ind_coupled(i, N, trial, 2)
            for test in 1:4
                test_ind_m = grid_to_ind_coupled(i, N, test, 1)
                test_ind_m1 = grid_to_ind_coupled(i, N, test , 2)
                for j in 1:gp

                    #going to assume v_A = 1?
                    #eq1 I
                    #I's are negative as I and W need to be opposite signs for eval problem
                    I[test_ind_m, trial_ind_m] -= -jac * Em[2, test, j] * r[j]^3 * Em[2, trial, j] * 1/jac^2 * wg[j] * dens[j]
                    I[test_ind_m, trial_ind_m] -= -jac * (m^2-1) * Em[1, test, j] * r[j] * Em[1, trial, j] * wg[j] * dens[j]
                    I[test_ind_m, trial_ind_m1] -= -jac * Em[2, test, j] * r[j]^3 * ϵ[j] * Em1[2, trial, j] * 1/jac^2 * wg[j] * dens[j]
                    #additional term due to derivative of density.
                    I[test_ind_m, trial_ind_m] -= jac * Em[1, test, j] * r[j]^2 * ddens[j] * Em[1, trial, j] * wg[j]

                    #damping contribution!
                    I[test_ind_m, trial_ind_m] -= jac * 1im * δ * (Em[3, test, j] / jac^2 + 1/r[j] * Em[2, test, j] / jac - m^2 / r[j]^2 * Em[1, test, j]) * dens[j] * (r[j] * Em[3, trial, j] / jac^2 + 3 * Em[2, trial, j] / jac + Em[1, trial, j] /r[j] - m^2/r[j] * Em[1, trial, j])

                    #eq1 W
                    W[test_ind_m, trial_ind_m] += jac * Em[2, test, j] * r[j]^3 * km[j]^2 * Em[2, trial, j] * 1/jac^2 * wg[j]
                    W[test_ind_m, trial_ind_m] += jac * (m^2-1) * Em[1, test, j] * r[j] * km[j]^2 * Em[1, trial, j] * wg[j]

                    #extra term that would be canceled by the current term.
                    #this implies that the current term is extremely important
                    #don't think we can just ignore it...
                    #yeah this completly changes everything!
                    #think re-creating continuum based of evals will help us bugfix our code though!
                    #W[test_ind_m, trial_ind_m] += jac * Em[1, test, j] * km[j] * (3*kmp[j] + r[j] * kmpp[j]) * Em[1, trial, j] * wg[j]
                    


                    #eq2 I
                    I[test_ind_m1, trial_ind_m1] -= -jac * Em1[2, test, j] * r[j]^3 * Em1[2, trial, j] * 1/jac^2 * wg[j] * dens[j]
                    I[test_ind_m1, trial_ind_m1] -= -jac * (m1^2-1) * Em1[1, test, j] * r[j] * Em1[1, trial, j] * wg[j] * dens[j]
                    I[test_ind_m1, trial_ind_m] -= -jac * Em1[2, test, j] * r[j]^3 * ϵ[j] * Em[2, trial, j] * 1/jac^2 * wg[j] * dens[j]
                    #additional term due to density derivative
                    #no negative is due to it failing with a negative, probably due to integration by parts!
                    I[test_ind_m1, trial_ind_m1] -= jac * Em1[1, test, j] * r[j]^2 * ddens[j] * Em1[1, trial, j] * wg[j]
                    
                    #damping cont.
                    I[test_ind_m1, trial_ind_m1] -= jac * 1im * δ * (Em1[3, test, j] / jac^2 + 1/r[j] * Em1[2, test, j] / jac - m1^2 / r[j]^2 * Em1[1, test, j]) * dens[j] * (r[j] * Em1[3, trial, j] / jac^2 + 3 * Em1[2, trial, j] / jac + Em1[1, trial, j] /r[j] - m1^2/r[j] * Em1[1, trial, j])

                    #eq2 W
                    W[test_ind_m1, trial_ind_m1] += jac * Em1[2, test, j] * r[j]^3 * km1[j]^2 * Em1[2, trial, j] * 1/jac^2 * wg[j]
                    W[test_ind_m1, trial_ind_m1] += jac * (m1^2-1) * Em1[1, test, j] * r[j] * km1[j]^2 * Em1[1, trial, j] * wg[j]

                    #extra term that would be canceled by the current term.
                    #W[test_ind_m1, trial_ind_m1] += jac * Em1[1, test, j] * km1[j] * (3*km1p[j] + r[j] * km1pp[j]) * Em1[1, trial, j] * wg[j]

                end
            end
        end
    end

    #apply the boundary conditions.
    if abs(m)==1
        bc1 = [2, 2*N-1]
    else
        bc1 = [1, 2*N-1]
    end
    if abs(m1)==1
        bc2 = [2*N+2, 4*N-1]
    else
        bc2 = [2*N+2, 4*N-1]
    end

    bcinds = vcat(bc1, bc2)

    for i in bcinds
        W[i, :] .= 0.0
        W[:, i] .= 0.0
        I[i, :] .= 0.0
        I[:, i] .= 0.0
    end

    for i in bcinds
        W[i, i] = 1.0
        I[i, i] = 1.0
    end 

    #even with current term, matrix seems `pretty` Hermitian, ie ~e-15 but now we get negative evals???
    #display(maximum(W - W')) #so these matrices are Hermitian to a tolerance of ~e-15
    #ω2, evals = eigen(Hermitian(W), Hermitian(I))
    ω2, evals = eigs(W, I, nev=5, ritzvec=true, sigma=(0.3259/R)^2)

    Emsol = evals[1:2:2*N, :]
    Em1sol = evals[2*N+1:2:end, :]
    return R .*sqrt.(ω2), Emsol, Em1sol

end

N = 4000
#rgrid = LinRange(0, 1, N);
#this makes a big difference!
rgrid = clustered_grid(N, 0.75, 0.8, 0.2)
ω, Em, Em1 = fudam(N, 10.0, 1, -1, rgrid, -4e-9);
#ω, ϕm, ϕm1 = phi_fudam(N, 10.0, 2, -2);

display(size(Em))

display(ω)

#2 from sqrt of R0, could be made up..
scatter(ones(length(ω)),  ω, ylimits=(-0.05, 1.05))

rdata = zeros(length(ω));
omdata = zeros(ComplexF64, length(ω));
col = zeros(length(ω));
for i in 1:1:length(ω)
    #this is much much better for cylinder case.
    #maybe this would be better in general then??
    rm = argmax(abs.(rgrid .*Em[:, i]))
    rm1 = argmax(abs.(rgrid .*Em1[:, i]))
    #rm = argmax(abs.(Em[:, i]))
    #rm1 = argmax(abs.(Em1[:, i]))
    if abs(Em[rm, i]) > abs(Em1[rm1, i])

        rdata[i] = rgrid[rm]
        col[i] = 0.0
    else
        rdata[i] = rgrid[rm1]
        col[i] = 1.0
    end

    omdata[i] = ω[i]
end

for i in 1:1:length(ω)
    #this is much much better for cylinder case.
    #maybe this would be better in general then??
    rm = argmax(abs.(ϕm[:, i]))
    rm1 = argmax(abs.(ϕm1[:, i]))
    #rm = argmax(abs.(Em[:, i]))
    #rm1 = argmax(abs.(Em1[:, i]))
    if abs(ϕm[rm, i]) > abs(ϕm1[rm1, i])

        rdata[i] = rgrid[rm]
        col[i] = 0.0
    else
        rdata[i] = rgrid[rm1]
        col[i] = 1.0
    end

    omdata[i] = ω[i]
end

#this kind of reconstructs the continuum, pretty garbage atm though!
#think the modes nearest the axis are just going to behave a little bit weird, may be an unavoidable problem.
#this is essentially caused by the boundary condition, these modes are somewhat artifically inflated at r=0
#scaling by r to get ϕ reduces the problem but gives a separate issue.
scatter(rdata, real.(omdata), ylimits=(-0.05, 1.05), group=col, markersize=4.0)#, legend=["m=1", "m=2"])

display(rdata[4])

plot(rgrid, Em[:, 25])

tae_ind = argmin(abs.(ω .- 0.319)) 
plot(rgrid, ϕm[:, tae_ind])
plot!(rgrid, ϕm1[:, tae_ind])
#little concerning that this is not exact..., hopefully this is due to the boundary conditions and not something else...
display(ω[tae_ind])

display(imag(ω[tae_ind])/real(ω[tae_ind]))
display(ω[1])

tae_ind = 1
#sacle is way wrong and eval is a bit small but otherwise pretty spot on!
#getting some stuff that resembles Bowden's case now we have density, but we have got multiple tae's in the frequency area. will need to improve resolution and/or add damping
#adding damping means we will probably need to change solver...
plot(rgrid, -real.(Em[:, tae_ind]))
plot!(rgrid, -real.(Em1[:, tae_ind]))
plot!(rgrid, -imag.(Em[:, tae_ind]))
plot!(rgrid, -imag.(Em1[:, tae_ind]))


plot(rgrid, rgrid .* Em[:, tae_ind])
plot!(rgrid, rgrid .* Em1[:, tae_ind])


plot(rgrid, rgrid .*Em[:, 800])

#lets see if we can plot a surface of this,

nθ = 50
θgrid = LinRange(0, 2π, nθ)


z = zeros(N, nθ, length(ω));

for i in 1:N

    for j in 1:nθ

        #z[i, j, :] = @. real( Em[i, :] * exp(1im * 1 * θgrid[j]) + Em1[i, :] * exp(1im * 2 * θgrid[j]))
        #think it is better to do this one to more easily connect with future work in terms of phi.
        #we may want to make sure we can solve this same equation but with phi instead of Em 
        #this will require checking some boundary conditions etc
        #also need to check the contribution of the current term.
        z[i, j, :] = @. real( rgrid[i] * Em[i, :] * exp(1im * 1 * θgrid[j]) + rgrid[i] * Em1[i, :] * exp(1im * 2 * θgrid[j]))
        #z[i, j, :] = @. real(ϕm[i, :] * exp(1im * 1 * θgrid[j]) + ϕm1[i, :] * exp(1im * 2 * θgrid[j]))
    end
end

scatter(rdata, omdata, ylimits=(-0.05, 1.05), group=col, markersize=4.0)

ind_m = argmin(abs.(ω .- 0.4795)) 
ind_m1 = argmin(abs.(ω .- 0.88))
tae_ind = argmin(abs.(ω .- 0.35)) 
#these two are either side of the tae.
combo_up_ind = argmin(abs.(ω .- 0.44)) 
combo_dn_ind = argmin(abs.(ω .- 0.28)) 

plot(rgrid, rgrid .* Em[:, ind_m])
plot!(rgrid, rgrid .*Em1[:, ind_m])

plot(rgrid, rgrid .* Em[:, ind_m1])
plot!(rgrid, rgrid .*Em1[:, ind_m1])

plot(rgrid, rgrid .* Em[:, tae_ind])
plot!(rgrid, rgrid .*Em1[:, tae_ind])

plot(rgrid, rgrid .* Em[:, combo_up_ind])
plot!(rgrid, rgrid .*Em1[:, combo_up_ind])

plot(rgrid, rgrid .* Em[:, combo_dn_ind])
plot!(rgrid, rgrid .*Em1[:, combo_dn_ind])



surface(θgrid, rgrid, z[:, :, ind_m])
surface(θgrid, rgrid, z[:, :, ind_m1])
surface(θgrid, rgrid, z[:, :, tae_ind])

#this mode is absolutly wild
#is essentially cos(x) - cos(2x) due to km = -km1
#but because of the singularity at r0, Em(r0-δ) = -Em(r0+δ) ie the function rapidly flips, so the total efunc is
#kind of like cos(x) - cos(2x) and -cos(x) + cos(2x) at the same time, but for slightly different r's so they 
#do not destroy each other! wild.
surface(θgrid, rgrid, z[:, :, combo_up_ind])

#this looks even more wild, for some reason for r<r0, only the cos(x) wave seems to be present
#for r>r0, we again get the sum of the two waves??
surface(θgrid, rgrid, z[:, :, combo_dn_ind])


#so for top one, both waves contribute equal magntiude on each side of r0 but with opposite sign, 
#for bottom one, only the m=1 wave has a meaningful contribution for r<r0, but for r>r0 both waves contribute equally and with the same sign?
#why does any of this happen??
#why do these flip sign over the singularity? I think this is very important because the `normal` waves far from the gap do not do this.


#also
#why are the perturbations confined to the flux surface but only kind of?
#why can the tae exist everywhere? Why is it not constrained??
#also why are there gaps lol

#in cyl limit, m=1 case is a constant until it hits the singularity, this seems weird
#need to understand physically why the m=1 mode is actually different at r=0 vs teh other modes.
#seems like they should be the same??


plot(rgrid, ϕm[:, tae_ind])
plot!(rgrid, ϕm1[:, tae_ind])

#try to plot ξ_r,

q = @. 1+ rgrid^2

#this have omega in it?? I guess that is just a scaling, we only really care about radial dependence.
ξ_r = @. (q * (-rgrid*(-1) + (1)*q)) * Em[:, ind_m] / (1+q^2)

plot(rgrid, ξ_r)
