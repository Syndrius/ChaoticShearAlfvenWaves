
"""
    label_mode(ϕft::Array{ComplexF64, 3}, grids::GridsT, rm::Array{Int64, 2}, ϕm::Array{Float64, 2})

Labels each eigenfunction with it's dominant mode number, the peak radial position of that mode.
"""
function label_mode(ϕft::Array{ComplexF64, 3}, grids::GridsT, rm::Array{Int64, 2}, ϕm::Array{Float64, 2})


    for j in 1:grids.θ.N, k in 1:grids.ζ.N

        rm[j, k] = argmax(abs.(real.(ϕft[:, j, k])))
        ϕm[j, k] = abs.(real.(ϕft[rm[j, k], j, k]))
    end

    max_mode = argmax(ϕm)

    mlab = mode_label(max_mode[1], grids.θ)
    nlab = mode_label(max_mode[2], grids.ζ)

    return rm[max_mode], (mlab, nlab)


end



"""
    label_mode(ϕft::Array{ComplexF64, 4}, grids::GridsT, rm::Array{Int64, 2}, ϕm::Array{Float64, 2})

Labels each eigenfunction with it's dominant mode number, the peak radial position of that mode.
"""
function label_mode(ϕft::Array{ComplexF64, 4}, grids::GridsT, rm::Array{Int64, 2}, ϕm::Array{Float64, 2})


    for j in 1:grids.θ.N, k in 1:grids.ζ.N

        rm[j, k] = argmax(abs.(real.(ϕft[:, j, k, 1])))
        ϕm[j, k] = abs.(real.(ϕft[rm[j, k], j, k, 1]))
    end

    max_mode = argmax(ϕm)

    mlab = mode_label(max_mode[1], grids.θ)
    nlab = mode_label(max_mode[2], grids.ζ)

    return rm[max_mode], (mlab, nlab)

end




#island case, subject to chaneg!!!!
function label_mode(ϕft::Array{ComplexF64, 3}, grids::MapGridsT, κm::Array{Int64, 2}, ϕm::Array{Float64, 2})

    #start from 2 to hopefully avoid the (0, 0) mode as a possibility??
    shift = 1 #to avoid (0, 0) modes.

    for j in 1+shift:grids.Nᾱ, k in 1:grids.Nφ

        κm[j-shift, k] = argmax(abs.(real.(ϕft[:, j, k])))
        ϕm[j-shift, k] = abs.(real.(ϕft[κm[j-shift, k], j, k]))
    end


    max_mode = argmax(ϕm)

    #uses different syntax for MapGridsT
    mlab, nlab = mode_label(max_mode[1] + shift, max_mode[2], grids)
    #mlab = mode_label(max_mode[1], grids.θ)
    #nlab = mode_label(max_mode[2], grids.ζ)

    return κm[max_mode], (mlab, nlab)

end