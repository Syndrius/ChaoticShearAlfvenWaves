
ϑgrid = MID.Structures.asm_grid(start=1, N=2)#, f_quad=1)
φgrid = MID.Structures.asm_grid(start=-1, N=1)#, f_quad=1)#, incr=2)
#grids = init_grids(N=N, mstart=1, mcount=2, nstart=-1, ncount=1);

#bounds chosen to not go outside the bounding surfaces
sgrid = MID.Structures.ContGridDataT(N=50, start=0.2, stop=0.8)
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
