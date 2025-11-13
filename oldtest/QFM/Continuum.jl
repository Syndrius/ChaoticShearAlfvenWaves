
ϑgrid = init_grid(type=:as, N=2, start=1)
φgrid = init_grid(type=:as, N=1, start=-1)
#grids = init_grids(N=N, mstart=1, mcount=2, nstart=-1, ncount=1);

#bounds chosen to not go outside the bounding surfaces
sgrid = init_grid(type=:rc, N=50, start=0.2, stop=0.8)
grids = init_grids(sgrid, ϑgrid, φgrid)


ω = compute_continuum(prob, grids, surfs);

gap = 0.32


mindiff = 10.0

for i in 1:ϑgrid.N, j in 1:φgrid.N
    diff = minimum(abs.(gap .- ω[:, i, j]))
    if diff < mindiff
        global mindiff = diff
    end
end

#this tests asserst that there is a gap in the continuum in the correct spot.
@test mindiff > 0.04
