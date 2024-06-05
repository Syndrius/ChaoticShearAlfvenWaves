

#direct implementation of Axel's equation for phi
#hopefully we can replicated Axel's results
#may not worry about doing this, seems a bit hard tbh!

#This is not even close! seems like a big problemo!

using FastGaussQuadrature
using MID
using LinearAlgebra
using Arpack
using SparseArrays
using Plots; plotlyjs()


function grid_to_ind_coupled(i, N, h, eqn)

    grid_id = [0, 0, 1, 1]

    basis_id = [0, 1, 0, 1]

    return 2 * (eqn-1)*N + 1 + 2*(i-1) + basis_id[h] + 2*grid_id[h]

end

#TODO!
function ϕm_to_cont(ω, N, ϕm, ϕm1)
    rdata = zeros(length(ω));
    omdata = zeros(ComplexF64, length(ω));
    col = zeros(length(ω));
    rgrid = clustered_grid(N, 0.93, 0.98, 0.25)
    for i in 1:1:length(ω)
        #this is much much better for cylinder case.
        #maybe this would be better in general then??
        rm = argmax(abs.(rgrid .*ϕm[:, i]))
        rm1 = argmax(abs.(rgrid .*ϕm1[:, i]))
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
    return rdata, omdata, col
end


function phi_2mode(N, m, R0, δ)
    #I think unfortunatly we have to do this now, think we are fked until we can replicated Axel's code.
    rgrid = clustered_grid(N, 0.93, 0.98, 0.25) #will just assume Axel's case for this.
    m1 = m+1
    n = -m
    gp = 5
    ξ, wg = gausslegendre(gp)
    H, dH, ddH = MID.Misc.hermite_basis(ξ)

    W = zeros(ComplexF64, 4*N, 4*N)
    I = zeros(ComplexF64, 4*N, 4*N)

    ϕm = zeros(3, 4, gp)
    ϕm1 = zeros(3, 4, gp)

    ϕm[1, :, :] = H
    ϕm[2, :, :] = dH
    ϕm[3, :, :] = ddH

    #probably not needed.
    ϕm1[1, :, :] = H
    ϕm1[2, :, :] = dH
    ϕm1[3, :, :] = ddH

    if abs(m)==0
        bc1 = [2, 2*N-1]
    else
        bc1 = [1, 2*N-1]
    end
    if abs(m1)==0
        bc2 = [2*N+2, 4*N-1]
    else
        bc2 = [2*N+1, 4*N-1]
    end

    bcinds = vcat(bc1, bc2)

    rows = Array{Int64}(undef, 0)
    cols = Array{Int64}(undef, 0)
    Wdata = Array{ComplexF64}(undef, 0)
    Idata = Array{ComplexF64}(undef, 0)

    for i in 1:N-1

        r, dr = MID.Misc.local_to_global(i, ξ, collect(rgrid))

        jac = dr/2

        Δp = r/(4*R0)

        Δpp = 1/(4*R0)

        ϵg = @. 2*Δp

        dϵg = 2 * Δpp * ones(length(r))

        ϵb = @. - r/R0 

        ϵ = @. ϵg/2 - 2 * ϵb # 5/2 * r/R good sign.

        dϵ = (Δpp - 2 / R0) * ones(length(r))

        q = @. 1.05 + 0.55*r^2

        dq = @. 2 * 0.55 * r

        dens = @. 1/2*(1-tanh((r-0.8)/0.1))

        ddens  = @. -5 * sech(8-10*r)^2

        #think we can do the weak form so that we don't need the derivative of km...
        km = @. (m/q + n)/R0 #will need to confirm this? may be -n instead
        km1 = @. (m1/q + n)/R0

        dk2m = @. -2 * m * (n + m/q) * dq / (R0^2*q^2)
        dk2m1 = @. -2 * m1 * (n + m1/q) * dq / (R0^2*q^2)

        dkmkm1 = @. -m1 * (n + m/q) * dq / (R0^2*q^2) -  m * (n + m1/q) * dq / (R0^2*q^2)


        for trial in 1:4

            trial_ind_m = grid_to_ind_coupled(i, N, trial, 1)
            trial_ind_m1 = grid_to_ind_coupled(i, N, trial, 2)
            for test in 1:4
                test_ind_m = grid_to_ind_coupled(i, N, test, 1)
                test_ind_m1 = grid_to_ind_coupled(i, N, test, 2)

                #I will be left hand side of Axel's equation.
                #weak form will defo multiply by r.
                #we will probably have to derive Axels result from equation given in paper.
                Isum11 = 0
                Isum12 = 0
                Isum21 = 0
                Isum22 = 0
                Wsum11 = 0
                Wsum12 = 0
                Wsum21 = 0
                Wsum22 = 0

                for j in 1:gp
                    
                    #ψ_m acting on equation 1.
                    ############################
                    #ψ_m ϕ_m terms
                    ######
                    #second deriv term
                    #-jac * ϕm[2, test, j] * r[j] * dens[j] * ϕm[2, trial, j] * wg[j] / jac^2
                    Isum11 += -jac * (ϕm[2, test, j] / r[j] / jac - ϕm[1, test, j] / r[j]^2) * r[j] * dens[j] * ϕm[2, trial, j] * wg[j] / jac

                    #m^2 term
                    Isum11 += -jac * ϕm[1, test, j] * m^2/r[j]^2 * dens[j] * ϕm[1, trial, j] * wg[j]

                    #second deriv term.
                    Wsum11 += -jac * (ϕm[2, test, j] / r[j] / jac - ϕm[1, test, j] / r[j]^2) * r[j] * km[j]^2 * ϕm[2, trial, j] * wg[j] / jac

                    #k2 deriv term
                    Wsum11 += -jac * ϕm[1, test, j] * dk2m[j] * ϕm[1, trial, j] * wg[j] / r[j]

                    #m^2 term
                    Wsum11 += -jac * ϕm[1, test, j] * m^2/r[j]^2 * km[j]^2 * ϕm[1, trial, j] * wg[j]

                    ##ψ_m ϕ_m1 terms
                    ######
                    #second deriv term
                    Isum12 += -jac * (ϕm[2, test, j] / r[j] / jac - ϕm[1, test, j] / r[j]^2) * r[j] * dens[j] * ϵ[j] * ϕm1[2, trial, j] * wg[j] / jac

                    #no deriv term
                    Isum12 += -jac * ϕm[1, test, j] * dens[j] * ϵ[j] * ϕm1[1, trial, j] * wg[j] / r[j]^2

                    #dens deriv term, split into two!
                    Isum12 += -jac * ϕm[1, test, j] * ddens[j] * ϵ[j] * ϕm1[1, trial, j] * wg[j] / r[j]
                    #eps deriv term.
                    Isum12 += -jac * ϕm[1, test, j] * dens[j] * dϵ[j] * ϕm1[1, trial, j] * wg[j] / r[j]

                    #second deriv term. Wrong sign perhap.
                    Wsum12 += -jac * (ϕm[2, test, j] / r[j] / jac - ϕm[1, test, j] / r[j]^2) * km[j] * km1[j] * ϵg[j]/2 * r[j] * ϕm1[2, trial, j] * wg[j] / jac

                    #no deriv term, Wrong sign perhap
                    Wsum12 += -jac * ϕm[1, test, j] * km[j] * km1[j] * ϵg[j]/2 * ϕm1[1, trial, j] * wg[j] / r[j]^2

                    #dens deriv term, split into two!
                    Wsum12 += -jac * ϕm[1, test, j] * dkmkm1[j] * ϵg[j]/2 * ϕm1[1, trial, j] * wg[j] / r[j]
                    #eps deriv.
                    Wsum12 += -jac * ϕm[1, test, j] * km[j] * km1[j] * dϵg[j]/2 * ϕm1[1, trial, j] * wg[j] / r[j]




                    #ψ_m1 acting on equation 2.
                    ############################
                    #ψ_m1 ϕ_m1 terms
                    ######
                    #second deriv term
                    Isum22 += -jac * (ϕm1[2, test, j] / r[j] / jac - ϕm1[1, test, j] / r[j]^2) * r[j] * dens[j] * ϕm1[2, trial, j] * wg[j] / jac

                    #m1^2 term
                    Isum22 += -jac * ϕm1[1, test, j] * m1^2/r[j]^2 * dens[j] * ϕm1[1, trial, j] * wg[j]

                    #second deriv term.
                    Wsum22 += -jac * (ϕm1[2, test, j] / r[j] / jac - ϕm1[1, test, j] / r[j]^2) * r[j] * km1[j]^2 * ϕm1[2, trial, j] * wg[j] / jac

                    #k2 deriv term
                    Wsum22 += -jac * ϕm1[1, test, j] * dk2m1[j] * ϕm1[1, trial, j] * wg[j] / r[j]

                    #m^2 term
                    Wsum22 += -jac * ϕm1[1, test, j] * m1^2/r[j]^2 * km1[j]^2 * ϕm1[1, trial, j] * wg[j]

                    ##ψ_m1 ϕ_m terms
                    ######
                    #second deriv term
                    Isum21 += -jac * (ϕm1[2, test, j] / r[j] / jac - ϕm1[1, test, j] / r[j]^2) * r[j] * dens[j] * ϵ[j] * ϕm[2, trial, j] * wg[j] / jac

                    #no deriv term
                    Isum21 += -jac * ϕm1[1, test, j] * dens[j] * ϵ[j] * ϕm[1, trial, j] * wg[j] / r[j]^2

                    #dens deriv term, split into two!
                    Isum21 += -jac * ϕm1[1, test, j] * ddens[j] * ϵ[j] * ϕm[1, trial, j] * wg[j] / r[j]
                    #eps deriv term.
                    Isum21 += -jac * ϕm1[1, test, j] * dens[j] * dϵ[j] * ϕm[1, trial, j] * wg[j] / r[j]

                    #second deriv term. Wrong sign perhap.
                    Wsum21 += -jac * (ϕm1[2, test, j] / r[j] / jac - ϕm1[1, test, j] / r[j]^2) * km[j] * km1[j] * ϵg[j]/2 * r[j] * ϕm[2, trial, j] * wg[j] / jac

                    #no deriv term, Wrong sign perhap
                    Wsum21 += -jac * ϕm1[1, test, j] * km[j] * km1[j] * ϵg[j]/2 * ϕm[1, trial, j] * wg[j] / r[j]^2

                    #dens deriv term, split into two!
                    #havent included jac for derivs of k etc yet, maybe they are needed!!!
                    Wsum21 += -jac * ϕm1[1, test, j] * dkmkm1[j] * ϵg[j]/2 * ϕm[1, trial, j] * wg[j] / r[j]
                    #eps deriv.
                    Wsum21 += -jac * ϕm1[1, test, j] * km[j] * km1[j] * dϵg[j]/2 * ϕm[1, trial, j] * wg[j] / r[j]

                    
                    #=
                    #ψ_m acting on equation 1.
                    ############################
                    #ψ_m ϕ_m terms
                    ######
                    #second deriv term
                    Isum11 += -jac * ϕm[2, test, j] * r[j] * dens[j] * ϕm[2, trial, j] * wg[j] / jac^2

                    #m^2 term
                    Isum11 += -jac * ϕm[1, test, j] * m^2/r[j] * dens[j] * ϕm[1, trial, j] * wg[j]

                    #second deriv term.
                    Wsum11 += -jac * ϕm[2, test, j] * r[j] * km[j]^2 * ϕm[2, trial, j] * wg[j] / jac^2

                    #k2 deriv term
                    Wsum11 += -jac * ϕm[1, test, j] * dk2m[j] * ϕm[1, trial, j] * wg[j]

                    #m^2 term
                    Wsum11 += -jac * ϕm[1, test, j] * m^2/r[j] * km[j]^2 * ϕm[1, trial, j] * wg[j]

                    ##ψ_m ϕ_m1 terms
                    ######
                    #second deriv term
                    Isum12 += -jac * ϕm[2, test, j] * r[j] * dens[j] * ϵ[j] * ϕm1[2, trial, j] * wg[j] / jac^2

                    #no deriv term
                    Isum12 += -jac * ϕm[1, test, j] * dens[j] * ϵ[j] * ϕm1[1, trial, j] * wg[j] / r[j]

                    #dens deriv term, split into two!
                    Isum12 += -jac * ϕm[1, test, j] * ddens[j] * ϵ[j] * ϕm1[1, trial, j] * wg[j] 
                    #eps deriv term.
                    Isum12 += -jac * ϕm[1, test, j] * dens[j] * dϵ[j] * ϕm1[1, trial, j] * wg[j] 

                    #second deriv term. Wrong sign perhap.
                    Wsum12 += -jac * ϕm[2, test, j] * km[j] * km1[j] * ϵg[j]/2 * r[j] * ϕm1[2, trial, j] * wg[j] / jac^2

                    #no deriv term, Wrong sign perhap
                    Wsum12 += -jac * ϕm[1, test, j] * km[j] * km1[j] * ϵg[j]/2 * ϕm1[1, trial, j] * wg[j] / r[j]

                    #dens deriv term, split into two!
                    Wsum12 += -jac * ϕm[1, test, j] * dkmkm1[j] * ϵg[j]/2 * ϕm1[1, trial, j] * wg[j]
                    #eps deriv.
                    Wsum12 += -jac * ϕm[1, test, j] * km[j] * km1[j] * dϵg[j]/2 * ϕm1[1, trial, j] * wg[j]




                    #ψ_m1 acting on equation 2.
                    ############################
                    #ψ_m1 ϕ_m1 terms
                    ######
                    #second deriv term
                    Isum22 += -jac * ϕm1[2, test, j] * r[j] * dens[j] * ϕm1[2, trial, j] * wg[j] / jac^2

                    #m1^2 term
                    Isum22 += -jac * ϕm1[1, test, j] * m1^2/r[j] * dens[j] * ϕm1[1, trial, j] * wg[j]

                    #second deriv term.
                    Wsum22 += -jac * ϕm1[2, test, j] * r[j] * km1[j]^2 * ϕm1[2, trial, j] * wg[j] / jac^2

                    #k2 deriv term
                    Wsum22 += -jac * ϕm1[1, test, j] * dk2m1[j] * ϕm1[1, trial, j] * wg[j]

                    #m^2 term
                    Wsum22 += -jac * ϕm1[1, test, j] * m1^2/r[j] * km1[j]^2 * ϕm1[1, trial, j] * wg[j]

                    ##ψ_m1 ϕ_m terms
                    ######
                    #second deriv term
                    Isum21 += -jac * ϕm1[2, test, j] * r[j] * dens[j] * ϵ[j] * ϕm[2, trial, j] * wg[j] / jac^2

                    #no deriv term
                    Isum21 += -jac * ϕm1[1, test, j] * dens[j] * ϵ[j] * ϕm[1, trial, j] * wg[j] / r[j]

                    #dens deriv term, split into two!
                    Isum21 += -jac * ϕm1[1, test, j] * ddens[j] * ϵ[j] * ϕm[1, trial, j] * wg[j] 
                    #eps deriv term.
                    Isum21 += -jac * ϕm1[1, test, j] * dens[j] * dϵ[j] * ϕm[1, trial, j] * wg[j] 

                    #second deriv term. Wrong sign perhap.
                    Wsum21 += -jac * ϕm1[2, test, j] * km[j] * km1[j] * ϵg[j]/2 * r[j] * ϕm[2, trial, j] * wg[j] / jac^2

                    #no deriv term, Wrong sign perhap
                    Wsum21 += -jac * ϕm1[1, test, j] * km[j] * km1[j] * ϵg[j]/2 * ϕm[1, trial, j] * wg[j] / r[j]

                    #dens deriv term, split into two!
                    #havent included jac for derivs of k etc yet, maybe they are needed!!!
                    Wsum21 += -jac * ϕm1[1, test, j] * dkmkm1[j] * ϵg[j]/2 * ϕm[1, trial, j] * wg[j]
                    #eps deriv.
                    Wsum21 += -jac * ϕm1[1, test, j] * km[j] * km1[j] * dϵg[j]/2 * ϕm[1, trial, j] * wg[j]
                    =#

                    ###########
                    #now consider damping terms
                    #damping term is added to the rhs, explains the signs.
                    #maybe we need to divide by r? to match other cases??
                    #start with the left double deriv terms
                    Isum11 += jac * 1im * δ * ϕm[3, test, j] * dens[j] * ϕm[3, trial, j] * wg[j] / jac^4
                    Isum11 += jac * 1im * δ * ϕm[3, test, j] * dens[j] * ϕm[2, trial, j] * wg[j] / r[j] / jac^3
                    Isum11 += jac * 1im * δ * ϕm[2, test, j] / jac * ( ddens[j] * m^2/r[j]^2 * ϕm[1, trial, j] - dens[j] * 2 * m^2/r[j]^3 * ϕm[1, trial, j] + dens[j] * m^2 / r[j]^2 * ϕm[2, trial, j] / jac) * wg[j]
                    #Isum11 += -jac * 1im * δ * ϕm[3, test, j] / jac^2 * dens[j]  * m^2/r[j]^2 * ϕm[1, test, j] * wg[j]

                    #now left single deriv terms
                    #fk me the weak form for this is cooked af.
                    Isum11 += jac * 1im * δ * (ϕm[2, test, j] / r[j] / jac - ϕm[1, test, j] / r[j]^2) * dens[j] * ϕm[3, trial, j] / jac^2 * wg[j]
                    Isum11 += jac * 1im * δ * (ϕm[2, test, j] / r[j] / jac - ϕm[1, test, j] / r[j]^2) * dens[j] * ϕm[2, trial, j] / r[j] / jac * wg[j]
                    Isum11 += -jac * 1im * δ * (ϕm[2, test, j] / r[j] / jac - ϕm[1, test, j] / r[j]^2) * dens[j] * ϕm[1, trial, j] * m^2 / r[j]^2 * wg[j]


                    #now no left deriv term, this has not been integrated by parts
                    Isum11 += -jac * 1im * δ * ϕm[1, test, j] * m^2/r[j]^2 * dens[j] * ϕm[3, trial, j] / jac^2 * wg[j]
                    Isum11 += -jac * 1im * δ * ϕm[1, test, j] * m^2/r[j]^2 * dens[j] * ϕm[2, trial, j] / jac / r[j] * wg[j]
                    Isum11 += jac * 1im * δ * ϕm[1, test, j] * m^4/r[j]^4 * dens[j] * ϕm[1, trial, j] * wg[j]


                    #now damping terms of second equation.
                    #start with the left double deriv terms
                    Isum22 += jac * 1im * δ * ϕm1[3, test, j] * dens[j] * ϕm1[3, trial, j] * wg[j] / jac^4 
                    Isum22 += jac * 1im * δ * ϕm1[3, test, j] * dens[j] * ϕm1[2, trial, j] * wg[j] / r[j] / jac^3 
                    Isum22 += jac * 1im * δ * ϕm1[2, test, j] / jac * ( ddens[j] * m1^2/r[j]^2 * ϕm1[1, trial, j] - dens[j] * 2 * m1^2/r[j]^3 * ϕm1[1, trial, j] + dens[j] * m1^2 / r[j]^2 * ϕm1[2, trial, j] / jac) * wg[j]

                    #now left single deriv terms
                    #fk me the weak form for this is cooked af.
                    Isum22 += jac * 1im * δ * (ϕm1[2, test, j] / r[j] / jac - ϕm1[1, test, j] / r[j]^2) * dens[j] * ϕm1[3, trial, j] / jac^2 * wg[j]
                    Isum22 += jac * 1im * δ * (ϕm1[2, test, j] / r[j] / jac - ϕm1[1, test, j] / r[j]^2) * dens[j] * ϕm1[2, trial, j] / r[j] / jac * wg[j]
                    Isum22 += -jac * 1im * δ * (ϕm1[2, test, j] / r[j] / jac - ϕm1[1, test, j] / r[j]^2) * dens[j] * ϕm1[1, trial, j] * m1^2 / r[j]^2 * wg[j]


                    #now no left deriv term, this has not been integrated by parts
                    Isum22 += -jac * 1im * δ * ϕm1[1, test, j] * m1^2/r[j]^2 * dens[j] * ϕm1[3, trial, j] / jac^2 * wg[j]
                    Isum22 += -jac * 1im * δ * ϕm1[1, test, j] * m1^2/r[j]^2 * dens[j] * ϕm1[2, trial, j] / jac / r[j] * wg[j]
                    Isum22 += jac * 1im * δ * ϕm1[1, test, j] * m1^4/r[j]^4 * dens[j] * ϕm1[1, trial, j] * wg[j]
                end

                
                push!(rows, test_ind_m)
                push!(cols, trial_ind_m)
                push!(Idata, Isum11)
                push!(Wdata, Wsum11)


                push!(rows, test_ind_m)
                push!(cols, trial_ind_m1)
                push!(Idata, Isum12)
                push!(Wdata, Wsum12)


                push!(rows, test_ind_m1)
                push!(cols, trial_ind_m)
                push!(Idata, Isum21)
                push!(Wdata, Wsum21)


                push!(rows, test_ind_m1)
                push!(cols, trial_ind_m1)
                push!(Idata, Isum22)
                push!(Wdata, Wsum22)

            end
        end
    end

    #still an awful way to do bc's but seems to work
    rowsbc = findall(x->x in bcinds, rows)
    colsbc = findall(x->x in bcinds, cols)
    Idata[rowsbc] .= 0.0
    Idata[colsbc] .= 0.0
    #rowsbcW = findall(x->x in bcinds, rowsW)
    #colsbcW = findall(x->x in bcinds, colsW)
    Wdata[rowsbc] .= 0.0
    Wdata[colsbc] .= 0.0
    for i in bcinds
        push!(rows, i)
        push!(cols, i)
        push!(Wdata, 1.0)
        push!(Idata, 1.0)
    end
    #this is fast as fk boi, now we have sparse fellas.
    W = sparse(rows, cols, Wdata)
    I = sparse(rows, cols, Idata)


    #then solve...
    #so looks like this form is not Hermitian, interesting, works if we don't specify Hermitian!
    #ω2, evals = eigen(Matrix(Hermitian(W)), Matrix(Hermitian(I)))
    #ω2, evals = eigen(Matrix(W), Matrix(I))
    #ω2, evals = eigs(W, I, nev=5, ritzvec=true, sigma=(0.389/R0)^2)
    #case for R0=20
    ω2, evals = eigs(W, I, nev=5, ritzvec=true, sigma=(0.395/R0)^2)

    ϕmsol = evals[1:2:2*N, :]
    ϕm1sol = evals[2*N+1:2:end, :]
    return R0 .*sqrt.(ω2), ϕmsol, ϕm1sol

        
end


#so W_12 and W_21 terms should be positive according to Axel's paper. With that assumption, we get
#0.3914205661164267 - 0.0008441257537601732im
#-0.0006608163608205249


#very unclear if Δ should be in this code or not. Paper deriving equation clearly doesn't have it
#but Axel's case does,  Above code has non-zero Δ so our result is probably a coincidennce :(
#our code with the cylindrical toroidal thing gives (without Δ!)
#0.39016064255026706 - 0.0008575887398034195im
#-0.0006691947475311519

#if we run the code with W_12, W_21 negative, which seems to match our derivation...
#we get, which seems to match the paper perfectly.
#0.3890088489360981 - 0.0011340161672794382im
#-0.0008822846478165998

#N = 1000, δ=7
#0.389020662787227 - 0.0012016681579530395im

#-0.0009349474865143953


#with R0=20, and signs to match paper results
#0.3949677161396188 - 0.0004987521042342625im
#-0.0003939819590584715

#using our code with test_met we get,
#0.39458550793479824 - 0.0005117760302768599im
#-0.0004038788097112989



#N=3000, δ=9, R0=10
#0.3890088489360979 - 0.0011340161672787237im
#-0.0008822846478160434
#R0=20
#0.394967716139619 - 0.0004987521042342667im
#-0.00039398195905847506

N = 3000
R0 = 20
ω_a, ϕm, ϕm1 = phi_2mode(N, 2, R0, -4.0e-9);


rdata, omdata, col = ϕm_to_cont(ω_a, N, ϕm, ϕm1);

#based on eye-test, actual tae frequency seems pretty spot on, but damping is way off.
#Flipping signs to match Axel's paper made the tae much closer to MID's prediction, but seem to be a bit further of the eye test of Axel's paper.
#either case has the damping completly wrong...
#looks like there is a typo in Axel's case, seems like we have replicated it, required changing the sign to match our derivation!
#next we need to see if we can derive Axels equation from ours, 
#then we should know what approximations are needed.
#damping is good now! yayzees.
scatter(rdata, real.(omdata).^2, ylimits=(-0.02, 0.25), group=col, markersize=4.0)

tae_ind = tae_ind = argmin(abs.(ω_a.^2 .- 0.15598)) 
tae_ind = 1

rgrid = clustered_grid(N, 0.93, 0.98, 0.25)
plot(rgrid, real.(ϕm[:, tae_ind]))
plot!(rgrid, real.(ϕm1[:, tae_ind]))

display(ω_a[tae_ind])
display(imag(ω_a[tae_ind]^2))

a=1 #stop us going to the bottom.




#terms before we included the extra / r
#= 
#ψ_m acting on equation 1.
                    ############################
                    #ψ_m ϕ_m terms
                    ######
                    #second deriv term
                    Isum11 += -jac * ϕm[2, test, j] * r[j] * dens[j] * ϕm[2, trial, j] * wg[j] / jac^2

                    #m^2 term
                    Isum11 += -jac * ϕm[1, test, j] * m^2/r[j] * dens[j] * ϕm[1, trial, j] * wg[j]

                    #second deriv term.
                    Wsum11 += -jac * ϕm[2, test, j] * r[j] * km[j]^2 * ϕm[2, trial, j] * wg[j] / jac^2

                    #k2 deriv term
                    Wsum11 += -jac * ϕm[1, test, j] * dk2m[j] * ϕm[1, trial, j] * wg[j]

                    #m^2 term
                    Wsum11 += -jac * ϕm[1, test, j] * m^2/r[j] * km[j]^2 * ϕm[1, trial, j] * wg[j]

                    ##ψ_m ϕ_m1 terms
                    ######
                    #second deriv term
                    Isum12 += -jac * ϕm[2, test, j] * r[j] * dens[j] * ϵ[j] * ϕm1[2, trial, j] * wg[j] / jac^2

                    #no deriv term
                    Isum12 += -jac * ϕm[1, test, j] * dens[j] * ϵ[j] * ϕm1[1, trial, j] * wg[j] / r[j]

                    #dens deriv term, split into two!
                    Isum12 += -jac * ϕm[1, test, j] * ddens[j] * ϵ[j] * ϕm1[1, trial, j] * wg[j] 
                    #eps deriv term.
                    Isum12 += -jac * ϕm[1, test, j] * dens[j] * dϵ[j] * ϕm1[1, trial, j] * wg[j] 

                    #second deriv term. Wrong sign perhap.
                    Wsum12 += -jac * ϕm[2, test, j] * km[j] * km1[j] * ϵg[j]/2 * r[j] * ϕm1[2, trial, j] * wg[j] / jac^2

                    #no deriv term, Wrong sign perhap
                    Wsum12 += -jac * ϕm[1, test, j] * km[j] * km1[j] * ϵg[j]/2 * ϕm1[1, trial, j] * wg[j] / r[j]

                    #dens deriv term, split into two!
                    Wsum12 += -jac * ϕm[1, test, j] * dkmkm1[j] * ϵg[j]/2 * ϕm1[1, trial, j] * wg[j]
                    #eps deriv.
                    Wsum12 += -jac * ϕm[1, test, j] * km[j] * km1[j] * dϵg[j]/2 * ϕm1[1, trial, j] * wg[j]




                    #ψ_m1 acting on equation 2.
                    ############################
                    #ψ_m1 ϕ_m1 terms
                    ######
                    #second deriv term
                    Isum22 += -jac * ϕm1[2, test, j] * r[j] * dens[j] * ϕm1[2, trial, j] * wg[j] / jac^2

                    #m1^2 term
                    Isum22 += -jac * ϕm1[1, test, j] * m1^2/r[j] * dens[j] * ϕm1[1, trial, j] * wg[j]

                    #second deriv term.
                    Wsum22 += -jac * ϕm1[2, test, j] * r[j] * km1[j]^2 * ϕm1[2, trial, j] * wg[j] / jac^2

                    #k2 deriv term
                    Wsum22 += -jac * ϕm1[1, test, j] * dk2m1[j] * ϕm1[1, trial, j] * wg[j]

                    #m^2 term
                    Wsum22 += -jac * ϕm1[1, test, j] * m1^2/r[j] * km1[j]^2 * ϕm1[1, trial, j] * wg[j]

                    ##ψ_m1 ϕ_m terms
                    ######
                    #second deriv term
                    Isum21 += -jac * ϕm1[2, test, j] * r[j] * dens[j] * ϵ[j] * ϕm[2, trial, j] * wg[j] / jac^2

                    #no deriv term
                    Isum21 += -jac * ϕm1[1, test, j] * dens[j] * ϵ[j] * ϕm[1, trial, j] * wg[j] / r[j]

                    #dens deriv term, split into two!
                    Isum21 += -jac * ϕm1[1, test, j] * ddens[j] * ϵ[j] * ϕm[1, trial, j] * wg[j] 
                    #eps deriv term.
                    Isum21 += -jac * ϕm1[1, test, j] * dens[j] * dϵ[j] * ϕm[1, trial, j] * wg[j] 

                    #second deriv term. Wrong sign perhap.
                    Wsum21 += -jac * ϕm1[2, test, j] * km[j] * km1[j] * ϵg[j]/2 * r[j] * ϕm[2, trial, j] * wg[j] / jac^2

                    #no deriv term, Wrong sign perhap
                    Wsum21 += -jac * ϕm1[1, test, j] * km[j] * km1[j] * ϵg[j]/2 * ϕm[1, trial, j] * wg[j] / r[j]

                    #dens deriv term, split into two!
                    #havent included jac for derivs of k etc yet, maybe they are needed!!!
                    Wsum21 += -jac * ϕm1[1, test, j] * dkmkm1[j] * ϵg[j]/2 * ϕm[1, trial, j] * wg[j]
                    #eps deriv.
                    Wsum21 += -jac * ϕm1[1, test, j] * km[j] * km1[j] * dϵg[j]/2 * ϕm[1, trial, j] * wg[j]
                =#