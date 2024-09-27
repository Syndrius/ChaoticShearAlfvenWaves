
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