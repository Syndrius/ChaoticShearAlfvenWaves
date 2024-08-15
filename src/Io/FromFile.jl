"""
    inputs_from_file(; dir::String)

Reads the prob and grids structs from file using JLD2.
"""
function inputs_from_file(; dir::String)

    grids = load_object(dir*"grids.jld2")
    prob = load_object(dir*"prob.jld2")
    
    return prob, grids

end


"""
    evals_from_file(; dir::String)

Reads the evals struct from file using JLD2.
"""
function evals_from_file(; dir::String)

    evals = load_object(dir*"evals.jld2")

    return evals
end


"""
    efunc_from_file(; dir::String, ind, ft=true)

Reads a single eigenfunction from file using JLD2. Set ft=false to read the non-ft eigenfunction.
"""
function efunc_from_file(; dir::String, ind, ft=true)

    if ft
        filename = dir * @sprintf("efuncs_ft/efunc%04d.jld2", ind)
    else
        filename = dir * @sprintf("efuncs/efunc%04d.jld2", ind)
    end

    ϕ = load_object(filename)

    return ϕ
end