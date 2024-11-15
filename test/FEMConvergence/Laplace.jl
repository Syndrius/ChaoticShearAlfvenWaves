
#not sure if this will be actually part of MID, or an adjacent part.
#here we more rigorously test convergence of our methods.

#comparing with 
#https://assets.cambridge.org/97805212/90494/excerpt/9780521290494_excerpt.pdf

#start with 2d laplace problemo.
#using Linear basis functions

using LinearAlgebra
using FastGaussQuadrature
using Plots


struct LB1d
    H :: Array{Float64, 2}
    dH :: Array{Float64, 2}
end

struct LB2d
    H :: Array{Float64, 4}
    dHx :: Array{Float64, 4}
    dHy :: Array{Float64, 4}
end

function compute_bcs(Nx, Ny)

    leftx = 1:Ny
    rightx = (Nx-1)*Ny:Nx*Ny

    

    lefty = 1:Ny:Nx*Ny


    righty = Nx:Ny:Nx*Ny

    #[x1y1, x1y2, ... x1, yN, x2y1, x2y2, .... xNy1, xNy2, ..., xNyN]

    return vcat(leftx, lefty, rightx, righty)

end

println(compute_bcs(4, 4))

function compute_bcs(N)

    

    return [1, N]

end


function grid_to_ind_linear_2d(xind::Int64, yind::Int64, hx, hy, Nx, Ny)

    xgrid = xind + hx - 1
    ygrid = yind + hy - 1

    #display(xgrid)

    if xgrid == Nx + 1
        #display("per")
        xgrid = 1
    end

    if ygrid == Ny + 1
        ygrid = 1
    end

    ind = Nx * (xgrid - 1) + (ygrid - 1) 

    return ind + 1

end


function grid_to_ind_linear_1d(xind::Int64, h, N)

    xgrid = xind + h - 1



    return xgrid 

end

function linear_basis(gp)

    S = LB1d(zeros(2, length(gp)), zeros(2, length(gp)))
    #don't think this is needed for linear.
    #t = @. (gp + 1)/(2) #converts to correct range for spline
    #H = zeros(2, length(gp))
    #dH = zeros(2, length(gp))

    S.H[1, :] = @. 1/2 * (1-gp)
    S.H[2, :] = @. 1/2 * (1+gp)

    S.dH[1, :] =  -1/2*ones(length(gp))
    S.dH[2, :] = 1/2*ones(length(gp))

    return S
end

function combine_basis(Sx, Sy, gpx, gpy)

    as = (2, 2, length(gpx), length(gpy))

    S = LB2d(zeros(as), zeros(as), zeros(as))

    for i in 1:length(gpx), j in 1:length(gpy)

        for y in 1:2, x in 1:2

            S.H[x, y, i, j] = Sx.H[x, i] * Sy.H[y, j]
            S.dHx[x, y, i, j] = Sx.dH[x, i] * Sy.H[y, j]
            S.dHy[x, y, i, j] = Sx.H[x, i] * Sy.dH[y, j]
        end
    end

    return S

end

function linear_basis(gpx, gpy)

    Sx = linear_basis(gpx)
    Sy = linear_basis(gpy)

    S = combine_basis(Sx, Sy, gpx, gpy)
    

    return S
end

function linear_local_to_global(xnode, ynode, ξx, ξy, xgrid, ygrid)
    
    if xnode == length(xgrid)
        dx = π + xgrid[1] - xgrid[xnode]
    else
        dx = xgrid[xnode+1] - xgrid[xnode]
    end

    if ynode == length(ygrid)
        dy = π + ygrid[1] - ygrid[ynode]
    else
        dy = ygrid[ynode+1] - ygrid[ynode]
    end

    mpx = @. (ξx + 1) /2 * dx
    mpy = @. (ξy + 1) /2 * dy

    xglobal = mpx .+ xgrid[xnode]
    yglobal = mpy .+ ygrid[ynode]

    return xglobal, yglobal, dx, dy
end


function linear_local_to_global(xnode, ξ, xgrid)
    

    dx = xgrid[xnode+1] - xgrid[xnode]
    


    mpx = @. (ξ + 1) /2 * dx

    xglobal = mpx .+ xgrid[xnode]

    return xglobal, dx
end


function laplace_1d(N::Int64)

    xgrid = LinRange(0, π, N)#+1)[1:end-1]


    gp = 5
    ξ, wg = gausslegendre(gp)

    S = linear_basis(ξ)

    M = zeros(Float64, N, N)

    A = zeros(Float64, N, N)

    bcs = compute_bcs(N)

    for i in 1:N-1


        x, dx = linear_local_to_global(i, ξ, xgrid)

        jac = dx / 2
        

        for trialx in 1:2

            r_ind = grid_to_ind_linear_1d(i, trialx, N)

            if r_ind in bcs
                continue
            end


            for testx in 1:2

                l_ind = grid_to_ind_linear_1d(i, testx, N)

                if l_ind in bcs
                    continue
                end


                for kx in 1:gp

                    M[l_ind, r_ind] += jac * (S.dH[testx, kx] * S.dH[trialx, kx]) / jac^2 * wg[kx] 

                    A[l_ind, r_ind] += jac * S.H[testx, kx] * S.H[trialx, kx] * wg[kx]

                end
            end
        end
    end

    for i in bcs
        M[i, i] = 1.0
        A[i, i] = 1.0
    end

    #display(M)
    #display(A)

    return eigen(Hermitian(M), Hermitian(A))
end

N = 30
evals, efuncs = laplace_1d(N);

println(evals)
println(sqrt.(evals))

xgrid = LinRange(0, π, N)#+1)[1:end-1];

plot(xgrid, efuncs[:, 28])

ind = argmax(abs.(efuncs[:, 28]))

sf = efuncs[ind, 28] / sin(26*xgrid[ind])

plot!(xgrid, sf*sin.(26 .* xgrid))

#t_evals = 1:1:N-2


plot(abs.(t_evals .- sqrt.(evals[3:end]) ))

function con_test(Nlist, k)

    #N_list = [4, 8, 16, 32, 64, 128, 256, 512]

    eval_list = zeros(length(Nlist))
    efunc_list = zeros(length(Nlist))



    for (i, N) in enumerate(N_list)



        evals, efuncs = laplace_1d(N);

        xgrid = LinRange(0, π, N)#+1)[1:end-1]

        #find a scale factor to make the amplitudes match

        efunc = efuncs[:, 2+k]

        ind = argmax(abs.(efunc))

        sf = efunc[ind] / sin(k*xgrid[ind])



        #may also fix the sign???
        truesol = @. sf * sin(k*xgrid)

        trueval = k

        #efunc_list[i] = sqrt.(sum((truesol .- efuncs[:, 2+k]) .^2))
        efunc_list[i] = maximum(abs.(truesol .- efuncs[:, 2+k]))

        #display(trueval)
        #display(sqrt(evals[2+k]))

        


        

        eval_list[i] = sqrt(evals[2+k]) - trueval
        """
        if efuncs[2, 2+k] > 0
            println(efuncs[:, 2+k])
            println(truesol)

            efunc_list[i] = sqrt.(sum((truesol .- efuncs[:, 2+k]) .^2))
        else
            println(efuncs[:, 2+k])
            println(-1 .* truesol)
            efunc_list[i] = sqrt.(sum((-1 .* truesol .- efuncs[:, 2+k]) .^2))
        end
        """
            
    end

    return eval_list, efunc_list

end

N_list = [8, 16, 32, 64, 128, 256]#, 512, 1024, 2048]

#N_list = [16]

ev, ef = con_test(N_list, 5);

#println(ef)

plot(N_list, ev)#, xaxis=:log10, yaxis=:log10)

#seems to converge to well...
#this doesn't really make sense, bit unfort.
#may need a more complicated problemo.
#unsure if increasing is just a floating point accuracy thing, very surprsied even N=8 
#is basically a perfect representation...
plot(N_list, ef)
display(ef[1])










function laplace_2d(Nx::Int64, Ny::Int64)

    #assume square NxN grid

    xgrid = LinRange(0, π, Nx)#+1)[1:end-1]
    ygrid = LinRange(0, π, Ny)#+1)[1:end-1]

    gpx = 5
    gpy = 5
    ξx, wgx = gausslegendre(gpx)
    ξy, wgy = gausslegendre(gpy)

    S = linear_basis(ξx, ξy)

    M = zeros(Float64, Nx*Ny, Nx*Ny)

    A = zeros(Float64, Nx*Ny, Nx*Ny)

    bcs = compute_bcs(Nx, Ny)

    

    for i in 1:Nx-1, j in 1:Ny-1

        x, y, dx, dy = linear_local_to_global(i, j, ξx, ξy, xgrid, ygrid)

        jac = dx * dy / 4

        

        for trialx in 1:2, trialy in 1:2

            r_ind = grid_to_ind_linear_2d(i, j, trialx, trialy, Nx, Ny)
            

            if r_ind in bcs
                continue
            end

            for testx in 1:2, testy in 1:2

                l_ind = grid_to_ind_linear_2d(i, j, testx, testy, Nx, Ny)

                if l_ind in bcs
                    continue
                end

                for kx in 1:gpx, ky in 1:gpy

                    M[l_ind, r_ind] += jac * (S.dHx[testx, testy, kx, ky] * S.dHx[trialx, trialy, kx, ky] * 2/dx * 2/dx + S.dHy[testx, testy, kx, ky] * S.dHy[trialx, trialy, kx, ky] * 2/dy * 2/dy) * wgx[kx] * wgy[ky]

                    #display(jac * (S.dHx[testx, testy, kx, ky] * S.dHx[trialx, trialy, kx, ky] + S.dHy[testx, testy, kx, ky] * S.dHy[trialx, trialy, kx, ky]) / jac^2 * wgx[kx] * wgy[ky])

                    A[l_ind, r_ind] += jac * S.H[testx, testy, kx, ky] * S.H[trialx, trialy, kx, ky] * wgx[kx] * wgy[ky]

                end
            end
        end
    end

    #display(S.H)
    for i in bcs
        M[i, i] = 1.0
        A[i, i] = 1.0
    end


    #println(M)


    evals, efuncs = eigen(M, A)

    bc_evs = length(unique(bcs))

    println(length(evals) - bc_evs)

    ef = zeros(length(evals)-bc_evs, Nx, Ny)



    for i in 1:length(evals) - bc_evs

        for j in 1:Nx, k in 1:Ny

            ef[i, j, k] = efuncs[grid_to_ind_linear_2d(j, k, 1, 1, Nx, Ny), i+bc_evs]
        end
    end

    #println(M)

    #println(bcs)

    return evals[1+bc_evs:end], ef

end

Nx = 8
Ny = Nx

evals, efuncs = laplace_2d(Nx, Ny);

#println(sqrt.(evals)[36:end])

xgrid = LinRange(0, π, Nx)#+1)[1:end-1]
ygrid = LinRange(0, π, Ny)#+1)[1:end-1] 


plot(ygrid, efuncs[2, 2, :])

contourf(ygrid, xgrid, efuncs[4, :, :])

println((evals))# ./ 20)


#plot(LinRange(0, π, length(efuncs[10, :])), efuncs[10, :])


function con_test_2d(Nlist, k, m, n)

    #will only really work with the non-degenerate cases, 
    #otherwise the efuncs get mixed together.
    #and will be very hard to split.

    #assume Ny = Nx

    #N_list = [4, 8, 16, 32, 64, 128, 256, 512]

    eval_list = zeros(length(Nlist))
    efunc_list = zeros(length(Nlist))



    for (i, N) in enumerate(Nlist)


        evals, efuncs = laplace_2d(N, N);

        xgrid = LinRange(0, π, N)#+1)[1:end-1]
        ygrid = LinRange(0, π, N)

        #find a scale factor to make the amplitudes match

        efunc = efuncs[k, :, :]

        ind = argmax(abs.(efunc))

        sf = efunc[ind] / (sin(m*xgrid[ind[1]])*sin(n*ygrid[ind[2]]))



        #may also fix the sign???
        truesol = zeros(N, N)
        for j in 1:N, l in 1:N
            truesol[j, l] = sf * sin(m*xgrid[j]) * sin(n*ygrid[l])
        end

        trueval = m^2 + n^2

        #efunc_list[i] = sqrt.(sum((truesol .- efuncs[:, 2+k]) .^2))
        efunc_list[i] = maximum(abs.(truesol .- efunc))

        #display(trueval)
        #display(sqrt(evals[2+k]))

        
        eval_list[i] = evals[k] - trueval

            
    end

    return eval_list, efunc_list

end

Nlist = [4, 8, 16, 32]

ev, ef = con_test_2d(Nlist, 4, 2, 2)

#maybe we need a more scientific method rather than just plotting!?

plot(Nlist, ev)#, xaxis=:log10, yaxis=:log10)

#same thing seems to be happening here
#v unclear, will need a better eval solver.
plot(Nlist, ef)

ξ, wg = gausslegendre(5)

println(wg)

S = linear_basis(ξ)
display(S.H[1, :])