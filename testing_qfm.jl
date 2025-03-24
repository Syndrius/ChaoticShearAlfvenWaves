
using Plots

#this is kinda slow lol.
#but looks to be pretty good.
surfs = construct_surfaces([5, 13, 8, 11, 3], [8, 21, 13, 18, 5])

#perfecto mundo.
plot_surfs(surfs)

surf_itp = create_interpolation(surfs)
#we will need to optimize this a fair bit.
#defs giving the exact same thing as Zhisongs version though
#just need to put all the peices toegther.
scos, tsin, ssin, tcos = action(13, 21);

#this is different to Zhisongs case, but similar enough to be satisfied
#think we probably need better interpolation
#or more surfaces
#or bounding surfaces
#or jacobia for solving
#but this should be adequate for continuing I think.
transform_coords(0.62, 2.0, 0.0, surf_itp)


#we are not interpolating the same thing.
display(surf_itp.scos_itp(0.62))

display(scos)
display(tsin)
display(ssin)
display(tcos)

Nθ = 100
θgrid = LinRange(0, 2*π, Nθ)
ζ = 0.0

α = zeros((Nθ, size(scos)[1], size(scos)[2]))

pqMpol = 24
pqNtor = 8
mlist = collect(range(0, pqMpol))

collect(-pqNtor:0)
nlist = [collect(0:pqNtor);collect(-pqNtor:-1)] 

for i in 1:Nθ
    for j in 1:size(scos)[1]
        for k in 1:size(scos)[2]
            α[i, j, k] = mlist[j] * θgrid[i] - nlist[k] * ζ
        end
    end
end


#α = @. 5 * θgrid - 8 * ζ



cosα = cos.(α)
sinα = sin.(α)

s = zeros(Nθ)
t = zeros(Nθ)

for i in 1:Nθ

    for j in 1:size(scos)[1]
        for k in 1:size(scos)[2]
            s[i] += scos[j, k] * cosα[i, j, k]
            t[i] += θgrid[i] + tsin[j, k] * sinα[i, j, k]
        end
    end
end



plot(t, s)