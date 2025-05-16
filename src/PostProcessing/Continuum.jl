
"""
    label_mode(ϕft::Array{ComplexF64, 3}, grids::GridsT, x1m::Array{Int64, 2}, ϕm::Array{Float64, 2})

Labels each eigenfunction with it's dominant mode number, the peak radial position of that mode.
"""
function label_mode(ϕft::Array{ComplexF64, 3}, grids::GridsT, x1m::Array{Int64, 2}, ϕm::Array{Float64, 2})


    for j in 1:grids.x2.N, k in 1:grids.x3.N

        x1m[j, k] = argmax(abs.(real.(ϕft[:, j, k])))
        ϕm[j, k] = abs.(real.(ϕft[x1m[j, k], j, k]))
    end

    max_mode = argmax(ϕm)

    mlab = mode_label(max_mode[1], grids.x2)
    nlab = mode_label(max_mode[2], grids.x3)

    return x1m[max_mode], (mlab, nlab)


end



"""
    label_mode(ϕft::Array{ComplexF64, 4}, grids::GridsT, x1m::Array{Int64, 2}, ϕm::Array{Float64, 2})

Labels each eigenfunction with it's dominant mode number, the peak radial position of that mode.
"""
function label_mode(ϕft::Array{ComplexF64, 4}, grids::GridsT, x1m::Array{Int64, 2}, ϕm::Array{Float64, 2})


    for j in 1:grids.x2.N, k in 1:grids.x3.N

        x1m[j, k] = argmax(abs.(real.(ϕft[:, j, k, 1])))
        ϕm[j, k] = abs.(real.(ϕft[x1m[j, k], j, k, 1]))
    end

    max_mode = argmax(ϕm)

    mlab = mode_label(max_mode[1], grids.x2)
    nlab = mode_label(max_mode[2], grids.x3)

    return x1m[max_mode], (mlab, nlab)

end

