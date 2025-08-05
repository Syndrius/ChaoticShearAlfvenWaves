using MID
using Plots
#then see if we can change the 2d cartesian case to only have derivs, so we can use our og code
#need to change to cylindrical coords and swap to periodic bc's.
#then we can test in 3d, 
#then test in parallel.
#%%
Nx = 25
Ny = 25

#seems to behave as expected.
#llam, lf = MID.Helmholtz.linear(Nx, Ny, 3, 3)
#%%

llam[15:20]

#with proper bc's the evals of this are now behaving as expected.
clam, cf = MID.Helmholtz.cubic(Nx, Ny, 4, 4);

println(clam[7])
clam
clam[20:30]
clam[4*Ny:4*Ny+5]
#%%
xgrid = LinRange(0, 1, Nx) 
ygrid = LinRange(0, 2π, Ny+1)[1:end-1]
clam[2]
#looks like this is ok with the extra tangent local to global jac.
#ideally we could actually test these results properly.
plot(xgrid, 2 .* real.(cf[24, :, 7, 1]))
plot!(xgrid, 2 .* real.(cf[13, :, 7, 3]))

plot(ygrid, real.(cf[2, 7, :, 1]))
plot!(ygrid, real.(cf[2, 7, :, 4]))
#%%
clam[8]
#not even kind of periodic...
plot(ygrid, real.(cf[4, 5, :]))
#ok so solutions are not really periodic and they look rubish.
contourf(ygrid, xgrid, real.(lf[4, :, :]))
contourf(ygrid, xgrid, real.(cf[2, :, :]))


#%%
Nlist = [8, 12, 16, 20, 24, 32]#, 64]#, 128]

ktest = [3, 4, 5, 6, 7, 8]
#ktrue = (ktest .- 2) .^ 2
ktrue = [5, 8, 10, 10, 13, 13]

lkerr = zeros(length(Nlist), length(ktest));
ckerr = zeros(length(Nlist), length(ktest));
lferr = zeros(length(Nlist), length(ktest));
cferr = zeros(length(Nlist), 3);
gp = 3
cfs = []
for (i, N) in enumerate(Nlist)
    xgrid = LinRange(0, π, N)
    llam, lf = MID.Helmholtz.linear(N, N, gp, gp)
    clam, cf = MID.Helmholtz.cubic(N, N, gp, gp)
    
    lkerr[i, :] = abs.(llam[ktest] .- ktrue)
    ckerr[i, :] = abs.(clam[ktest] .- ktrue)

    truef = true_f(N, N);
    display(size(truef))
    display(size(cf))
    sf1 = truef[1, 3, 3] / cf[1, 3, 3]
    sf2 = truef[2, 3, 3] / cf[4, 3, 3]
    sf3 = truef[3, 3, 3] / cf[11, 3, 3]
    cferr[i, 1] = mean(truef[1, :, :] .- sf1 .* cf[1, :, :])
    cferr[i, 2] = mean(truef[2, :, :] .- sf2 .* cf[4, :, :])
    cferr[i, 3] = mean(truef[3, :, :] .- sf3 .* cf[11, :, :])
end
#%%
#cool, eigensolutions look good!
p = plot()
for k in 1:3
    scatter!(cferr[:, k])
end
display(p)
#%%
p = plot()#ylimits=(0, 1))
for i in 1:length(ktest)
    #scatter!((lkerr[:, i]))
    #scatter!((ckerr[:, i]))
    #scatter!(log.(lkerr[:, i]))
    #scatter!(log.(ckerr[:, i]))
    #scatter!((lferr[2:end, i]))
    #scatter!((cferr[1:end, i]))
end

display(p)
#%%

function true_f(Nx, Ny)
    xgrid = LinRange(0, π, Nx)
    ygrid = LinRange(0, π, Ny)
    truef = zeros(3, Nx, Ny)

    for i in 1:Nx
        for j in 1:Ny
            truef[1, i, j] = sin(1*xgrid[i]) * sin(1*ygrid[j])
            truef[2, i, j] = sin(2*xgrid[i]) * sin(2*ygrid[j])
            truef[3, i, j] = sin(3*xgrid[i]) * sin(3*ygrid[j])
            #truef[2, i, j] = sin(1*xgrid[i]) * sin(2*ygrid[j])
            #truef[3, i, j] = sin(2*xgrid[i]) * sin(1*ygrid[j])
            #truef[4, i, j] = sin(2*xgrid[i]) * sin(2*ygrid[j])
            #truef[5, i, j] = sin(1*xgrid[i]) * sin(3*ygrid[j])
            #truef[6, i, j] = sin(3*xgrid[i]) * sin(1*ygrid[j])
            #truef[7, i, j] = sin(2*xgrid[i]) * sin(3*ygrid[j])
        end
    end

    return truef
end
#%%
#for error, we can probably only really check the symmetric ones unless we make some complex condiftions.



#%%
Nx = 15
Ny = 15
xgrid = LinRange(0, π, Nx)
ygrid = LinRange(0, π, Ny)

truef = true_f(Nx, Ny);
#%%
contourf(ygrid, xgrid, truef[3, :, :])
contourf(ygrid, xgrid, cf[11, :, :])


#%%
truef = sin.(1 .*xgrid) .* sin.(2 .* ygrid')
plot(xgrid, lf[2, :, 3])
size(lf)
plot(xgrid, lf[4:Ny:end, 3])

plot(xgrid, truef[4, :])
contourf(ygrid, xgrid, truef[5, :, :])
#linear index to grid is probably cooked. oh well.
contourf(ygrid, xgrid, lf[5, :, :])
size(cf)
contourf(ygrid, xgrid, cf[3, :, :])

plot(xgrid, lf[4, :, 4])
llam

