
#damping just solving for Em for some clarity.


using Plots; plotlyjs()
using FastGaussQuadrature
using MID
using LinearAlgebra
using Arpack
using SparseArrays


function grid_to_ind_coupled(i, N, h, eqn)

    grid_id = [0, 0, 1, 1]

    basis_id = [0, 1, 0, 1]

    return 2 * (eqn-1)*N + 1 + 2*(i-1) + basis_id[h] + 2*grid_id[h]

end


function Em_to_cont(ω, rgrid, Em, Em1)
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
    return rdata, omdata, col
end

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

    if abs(m)==1
        bc1 = [2, 2*N-1]
    else
        bc1 = [1, 2*N-1]
    end
    if abs(m1)==1
        bc2 = [2*N+2, 4*N-1]
    else
        bc2 = [2*N+1, 4*N-1]
    end

    bcinds = vcat(bc1, bc2)
    #display(bcinds)
    rowsW = Array{Int64}(undef, 0)
    colsW = Array{Int64}(undef, 0)
    Wdata = Array{ComplexF64}(undef, 0)
    Idata = Array{ComplexF64}(undef, 0)
    rowsI = Array{Int64}(undef, 0)
    colsI = Array{Int64}(undef, 0)
    

    for i in 1:N-1

        r, dr = MID.Misc.local_to_global(i, ξ, collect(rgrid))

        #Berk's case,
        #=
        q = @. 1 + r^2

        dens = ones(length(r))
        ddens = zeros(length(r))
        =#
        

        #Bowden singular
        #=
        q = @. 1.0 + 2.0*r^2

        dens = @. 1/2*(1-tanh((r-0.7)/0.05))

        ddens = @. -10 * sech((r-0.07)/0.05)^2
        
        =#
        #Axel 2012
        
        q = @. 1.05 + 0.55*r^2

        #dens = ones(length(r)) #for testing.
        #ddens = zeros(length(r))
        dens = @. 1/2*(1-tanh((r-0.8)/0.1))
        ddens = @. -5 * sech(8-10*r)^2
        

        qpp = 2*ones(length(r))

        km = @. 1/R*(m/q + n)

        #kmp = @. -m*qp/(R*q^2)

        #kmpp = @. m * (2*qp^2 - q*qpp)/(R*q^3)

        km1 = @. 1/R*(m1/q + n)

        #km1p = @. -m1*qp/(R*q^2)

        #km1pp = @. m1 * (2*qp^2 - q*qpp)/(R*q^3)

        jac = dr / 2

        ϵ = @. 5*r/(2*R)
        #Fu Dam case 
        #this gives the exact frequency that Berk quotes in fig 1...
        #ϵ = @. 3*r/(2*R)

        #cylindrical limit!
        #ϵ = zeros(length(r))

        


        #equation 1
        for trial in 1:4

            trial_ind_m = grid_to_ind_coupled(i, N, trial, 1)
            trial_ind_m1 = grid_to_ind_coupled(i, N, trial, 2)
            #=
            if trial_ind_m in bcinds
                continue
            end
            if trial_ind_m1 in bcinds
                continue
            end
            =#
            for test in 1:4
                test_ind_m = grid_to_ind_coupled(i, N, test, 1)
                test_ind_m1 = grid_to_ind_coupled(i, N, test , 2)
                #=
                if test_ind_m in bcinds
                    continue
                end
                if test_ind_m1 in bcinds
                    continue
                end
                =#
                for j in 1:gp

                    push!(rowsI, test_ind_m)
                    push!(colsI, trial_ind_m)
                    push!(Idata, jac * Em[2, test, j] * r[j]^3 * Em[2, trial, j] * 1/jac^2 * wg[j] * dens[j])

                    push!(rowsI, test_ind_m)
                    push!(colsI, trial_ind_m)
                    push!(Idata, jac * (m^2-1) * Em[1, test, j] * r[j] * Em[1, trial, j] * wg[j] * dens[j])

                    push!(rowsI, test_ind_m)
                    push!(colsI, trial_ind_m)
                    push!(Idata, -jac * Em[1, test, j] * r[j]^2 * ddens[j] * Em[1, trial, j] * wg[j])# * jac)

                    push!(rowsI, test_ind_m)
                    push!(colsI, trial_ind_m1)
                    push!(Idata, jac * Em[2, test, j] * r[j]^3 * ϵ[j] * Em1[2, trial, j] * 1/jac^2 * wg[j] * dens[j])

                    push!(rowsI, test_ind_m)
                    push!(colsI, trial_ind_m)
                    push!(Idata, -jac * 1im * δ * (Em[3, test, j] / jac^2 + 1/r[j] * Em[2, test, j] / jac - m^2 / r[j]^2 * Em[1, test, j]) * dens[j] * (r[j] * Em[3, trial, j] / jac^2 + 3 * Em[2, trial, j] / jac + Em[1, trial, j] /r[j] - m^2/r[j] * Em[1, trial, j]))

                    #going to assume v_A = 1?
                    #eq1 I
                    #I's are negative as I and W need to be opposite signs for eval problem
                    #=
                    I[test_ind_m, trial_ind_m] -= -jac * Em[2, test, j] * r[j]^3 * Em[2, trial, j] * 1/jac^2 * wg[j] * dens[j]
                    I[test_ind_m, trial_ind_m] -= -jac * (m^2-1) * Em[1, test, j] * r[j] * Em[1, trial, j] * wg[j] * dens[j]
                    I[test_ind_m, trial_ind_m1] -= -jac * Em[2, test, j] * r[j]^3 * ϵ[j] * Em1[2, trial, j] * 1/jac^2 * wg[j] * dens[j]
                    #additional term due to derivative of density.
                    I[test_ind_m, trial_ind_m] -= jac * Em[1, test, j] * r[j]^2 * ddens[j] * Em[1, trial, j] * wg[j]

                    #damping contribution!
                    I[test_ind_m, trial_ind_m] -= jac * 1im * δ * (Em[3, test, j] / jac^2 + 1/r[j] * Em[2, test, j] / jac - m^2 / r[j]^2 * Em[1, test, j]) * dens[j] * (r[j] * Em[3, trial, j] / jac^2 + 3 * Em[2, trial, j] / jac + Em[1, trial, j] /r[j] - m^2/r[j] * Em[1, trial, j])
                    =#
                    push!(rowsW, test_ind_m)
                    push!(colsW, trial_ind_m)
                    push!(Wdata, jac * Em[2, test, j] * r[j]^3 * km[j]^2 * Em[2, trial, j] * 1/jac^2 * wg[j])

                    push!(rowsW, test_ind_m)
                    push!(colsW, trial_ind_m)
                    push!(Wdata, jac * (m^2-1) * Em[1, test, j] * r[j] * km[j]^2 * Em[1, trial, j] * wg[j]) 
                    #eq1 W
                    #=
                    W[test_ind_m, trial_ind_m] += jac * Em[2, test, j] * r[j]^3 * km[j]^2 * Em[2, trial, j] * 1/jac^2 * wg[j]
                    W[test_ind_m, trial_ind_m] += jac * (m^2-1) * Em[1, test, j] * r[j] * km[j]^2 * Em[1, trial, j] * wg[j]
                    =#

                    #extra term that would be canceled by the current term.
                    #this implies that the current term is extremely important
                    #don't think we can just ignore it...
                    #yeah this completly changes everything!
                    #think re-creating continuum based of evals will help us bugfix our code though!
                    #W[test_ind_m, trial_ind_m] += jac * Em[1, test, j] * km[j] * (3*kmp[j] + r[j] * kmpp[j]) * Em[1, trial, j] * wg[j]
                    push!(rowsI, test_ind_m1)
                    push!(colsI, trial_ind_m1)
                    push!(Idata, jac * Em1[2, test, j] * r[j]^3 * Em1[2, trial, j] * 1/jac^2 * wg[j] * dens[j])

                    push!(rowsI, test_ind_m1)
                    push!(colsI, trial_ind_m1)
                    push!(Idata, jac * (m1^2-1) * Em1[1, test, j] * r[j] * Em1[1, trial, j] * wg[j] * dens[j])

                    push!(rowsI, test_ind_m1)
                    push!(colsI, trial_ind_m1)
                    push!(Idata, -jac * Em1[1, test, j] * r[j]^2 * ddens[j] * Em1[1, trial, j] * wg[j])# * jac)

                    push!(rowsI, test_ind_m1)
                    push!(colsI, trial_ind_m)
                    push!(Idata, jac * Em1[2, test, j] * r[j]^3 * ϵ[j] * Em[2, trial, j] * 1/jac^2 * wg[j] * dens[j])

                    push!(rowsI, test_ind_m1)
                    push!(colsI, trial_ind_m1)
                    push!(Idata, -jac * 1im * δ * (Em1[3, test, j] / jac^2 + 1/r[j] * Em1[2, test, j] / jac - m1^2 / r[j]^2 * Em1[1, test, j]) * dens[j] * (r[j] * Em1[3, trial, j] / jac^2 + 3 * Em1[2, trial, j] / jac + Em1[1, trial, j] /r[j] - m1^2/r[j] * Em1[1, trial, j]))


                    #eq2 I
                    #=
                    I[test_ind_m1, trial_ind_m1] -= -jac * Em1[2, test, j] * r[j]^3 * Em1[2, trial, j] * 1/jac^2 * wg[j] * dens[j]
                    I[test_ind_m1, trial_ind_m1] -= -jac * (m1^2-1) * Em1[1, test, j] * r[j] * Em1[1, trial, j] * wg[j] * dens[j]
                    I[test_ind_m1, trial_ind_m] -= -jac * Em1[2, test, j] * r[j]^3 * ϵ[j] * Em[2, trial, j] * 1/jac^2 * wg[j] * dens[j]
                    #additional term due to density derivative
                    #no negative is due to it failing with a negative, probably due to integration by parts!
                    I[test_ind_m1, trial_ind_m1] -= jac * Em1[1, test, j] * r[j]^2 * ddens[j] * Em1[1, trial, j] * wg[j]
                    
                    #damping cont.
                    I[test_ind_m1, trial_ind_m1] -= jac * 1im * δ * (Em1[3, test, j] / jac^2 + 1/r[j] * Em1[2, test, j] / jac - m1^2 / r[j]^2 * Em1[1, test, j]) * dens[j] * (r[j] * Em1[3, trial, j] / jac^2 + 3 * Em1[2, trial, j] / jac + Em1[1, trial, j] /r[j] - m1^2/r[j] * Em1[1, trial, j])
                    =#

                    push!(rowsW, test_ind_m1)
                    push!(colsW, trial_ind_m1)
                    push!(Wdata, jac * Em1[2, test, j] * r[j]^3 * km1[j]^2 * Em1[2, trial, j] * 1/jac^2 * wg[j])

                    push!(rowsW, test_ind_m1)
                    push!(colsW, trial_ind_m1)
                    push!(Wdata, jac * (m1^2-1) * Em1[1, test, j] * r[j] * km1[j]^2 * Em1[1, trial, j] * wg[j]) 

                    #eq2 W
                    #=
                    W[test_ind_m1, trial_ind_m1] += jac * Em1[2, test, j] * r[j]^3 * km1[j]^2 * Em1[2, trial, j] * 1/jac^2 * wg[j]
                    W[test_ind_m1, trial_ind_m1] += jac * (m1^2-1) * Em1[1, test, j] * r[j] * km1[j]^2 * Em1[1, trial, j] * wg[j]
                    =#

                    #extra term that would be canceled by the current term.
                    #W[test_ind_m1, trial_ind_m1] += jac * Em1[1, test, j] * km1[j] * (3*km1p[j] + r[j] * km1pp[j]) * Em1[1, trial, j] * wg[j]

                end
            end
        end
    end

    #apply the boundary conditions.
    #unsure if this needs to be set or not, based on bowden singular result, it shouldn't be...
    #=
    for i in bcinds
        push!(rowsI, i)
        push!(colsI, i)
        push!(rowsW, i)
        push!(colsW, i)
        push!(Wdata, 1.0)
        push!(Idata, 1.0)
    end
    =#
    #this may be slow af, but should work as an awful solution...
    #hack solution to keep the matrices sparse, works, but probably not very well.
    rowsbcI = findall(x->x in bcinds, rowsI)
    colsbcI = findall(x->x in bcinds, colsI)
    Idata[rowsbcI] .= 0.0
    Idata[colsbcI] .= 0.0
    rowsbcW = findall(x->x in bcinds, rowsW)
    colsbcW = findall(x->x in bcinds, colsW)
    Wdata[rowsbcW] .= 0.0
    Wdata[colsbcW] .= 0.0
    for i in bcinds
        push!(rowsI, i)
        push!(colsI, i)
        push!(rowsW, i)
        push!(colsW, i)
        push!(Wdata, 1.0)
        push!(Idata, 1.0)
    end
    #this is fast as fk boi, now we have sparse fellas.
    W = sparse(rowsW, colsW, Wdata)
    I = sparse(rowsI, colsI, Idata)
    #=
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
    =#
    

    #display(Matrix(I))
    #display(I)
    #even with current term, matrix seems `pretty` Hermitian, ie ~e-15 but now we get negative evals???
    #display(maximum(W - W')) #so these matrices are Hermitian to a tolerance of ~e-15
    ω2, evals = eigen(Matrix(Hermitian(W)), Matrix(Hermitian(I)))
    #ω2, evals = eigs(W, I, nev=5, ritzvec=true, sigma=(0.332/R)^2)
    #ω2, evals = eigs(W, I, nev=5, ritzvec=true, sigma=(0.39/R)^2)

    Emsol = evals[1:2:2*N, :]
    Em1sol = evals[2*N+1:2:end, :]
    return R .*sqrt.(ω2), Emsol, Em1sol

end


#so 3000 and -4e-8 gives -8.64e-4, pretty close to Axel's result. 
#we seem to be converging to a slightly lower damping rate, this could probably be attributed to 
#less interaction of the neighbouring modes.
#we may have to expand this to consider the next two neighbouring modes...

#this example may show that our code is actually working!
#For R0=10, this produces ~-8.3e-4, (Ballpark is Axel's result but a bit lower) and MID produces ~-1.5e-3, almost double the damping. (note that this was done with 2 modes!) But when we move to R0=20, so that additional toroidal effects that we are included in MID are reduced, this prodices ~-5.4e-4, while MID produces ~-5.9e-4, so the gap is significanlty reduced.

#obviously further verification is required. But this is positive!
#we probably need to consider cases with more modes, which will require re-writing this...
#ideally we would get the bowden singular one to work,
#but density seems to cook it completly.
#with extra modes in MID, ie m=1 and m=4, we get ~-5.76e-4 so even closer together.

N = 100
#rgrid = LinRange(0, 1, N);
#this makes a big difference!
rgrid = clustered_grid(N, 0.92, 0.98, 0.2)
#rgrid = clustered_grid(N, 0.75, 0.85, 0.3)
ω, Em, Em1 = fudam(N, 10.0, 2, -2, rgrid, -0.0e-7);


rdata, omdata, col = Em_to_cont(ω, rgrid, Em, Em1);


scatter(rdata, real.(omdata), ylimits=(-0.05, 1.05), group=col, markersize=4.0)

tae_ind = argmin(abs.(omdata .- 0.387)) 
tae_ind = 1
#something is still wrong with the damping!!
plot(rgrid, -real.(Em[:, tae_ind]))
plot!(rgrid, -real.(Em1[:, tae_ind]))
plot!(rgrid, -imag.(Em[:, tae_ind]))
plot!(rgrid, -imag.(Em1[:, tae_ind]))

display(ω[tae_ind])
display(imag(ω[tae_ind])/real(ω[tae_ind]))
display(imag(ω[tae_ind]^2))


#fudam result is ω = 0.93 * (km * v_A)_{q=1.5}

#our km is defined differently,
#km = (m/q+n)/R
#gap appears at
#km = -km1
#(m/q+n) = ((m+1)/q+n) -> 1/q - 1 = -2/q + 1
#-> 3/q = 2 -> q=3/2 = 1.5, same.
#km_{q=1.5} = 2/3-1 = -1/3
#v_A{q=1.5} = surely has to be at the axis? otherwise this will depend on θ?
#then v_A = 1?
#so fu_dam result should be 0.93/3 = 0.31. This is what Berk states...
#think Berk may juts be using FuDams result???
#pretty sure 0.31 from q=1+r^2 is a fudam result, so no Δ


using Printf
using DelimitedFiles
using Plots; plotlyjs()

#would be quite nice to have an automatic convergence test function...
#ie just pass in a couple of lists and it does it all??
#considering we will be doing that often?
Nlist = [200, 500, 1000, 1500, 2000, 3000, 4000]
dlist = [-4e-7, -4e-8, -4e-9, -4e-10, -4e-11, -4e-12]
dlabs = [7, 8, 9, 10, 11, 12]


dir_base = "data/fu_dam_singular_convergence_20/"

geo = GeoParamsT(R0=20.0)
#this took a while... probably should be using sparse matrices.
#axel case converges to ~8.19e-4, so underestimating the damping...
#could be the difference between assuming km=-km1?
#may need to check if these two get closer together as R0 gets larger.

#Bowden singular case
#converges to ~0.0174, compared to quoted result from singular method of 0.0176, and quoted method of
#contour method wich converges to ~0.0175, so pretty spot on.
for N in Nlist
    rgrid = clustered_grid(N, 0.75, 0.8, 0.2)
    #grids = init_grids(rgrid=rgrid, mstart=1, mcount=4, nstart=-2, ncount=1);
    for (i, δ) in enumerate(dlist)
        #prob = init_problem(q=Axel_q, geo=geo, δ=δ, dens=axel_dens);
        #prob = init_problem(q=singular_bowden_q, geo=geo, δ=δ, dens=bowden_singular_dens);

        dir = dir_base * @sprintf("N_%s_delta_%s/", N, dlabs[i])
        mkdir(dir)

        #inputs_to_file(prob=prob, grids=grids, dir=dir)

        #ω, ϕ = construct_and_solve(prob=prob, grids=grids, full_spectrum=false, σ=(0.391/geo.R0)^2, reconstruct=false);
        
        ω, Em, Em1 = fudam(N, 20.0, 1, -1, rgrid, δ);

        #just going to have to assume ind=1 is correct.
        eigvals_to_file(ω=ω, filename=dir*"vals.dat")

        #eigfuncs_to_file(ϕ=ϕ, filename=dir*"funcs.dat")

    end
end
