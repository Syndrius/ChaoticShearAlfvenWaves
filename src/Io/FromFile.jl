"""
    inputs_from_file(; dir::String)

Reads the prob and grids structs from file using JLD2.
"""
function inputs_from_file(dir::String)

    grids = load_object(joinpath(dir,"grids.jld2"))
    prob = load_object(joinpath(dir,"prob.jld2"))
    solver = load_object(joinpath(dir,"solver.jld2"))
    
    return prob, grids, solver

end


"""
    evals_from_file(; dir::String)

Reads the evals struct from file using JLD2.
"""
function evals_from_file(dir::String)

    evals = load_object(joinpath(dir,"evals.jld2"))

    return evals
end


"""
    efunc_from_file(; dir::String, ind::Int64, ft::Bool=true, deriv::Bool=false)

Reads a single eigenfunction from file using JLD2. Set ft=false to read the non-ft eigenfunction.
"""
function efunc_from_file(dir::String, ind::Int64; ft::Bool=true, deriv::Bool=false)

    if ft
        if deriv
            filename = joinpath(dir, @sprintf("efuncs_ft_deriv/efunc%05d.jld2", ind))
        else
            filename = joinpath(dir, @sprintf("efuncs_ft/efunc%05d.jld2", ind))
        end
    elseif deriv
        filename = joinpath(dir, @sprintf("efuncs_deriv/efunc%05d.jld2", ind))
    else
        filename = joinpath(dir, @sprintf("efuncs/efunc%05d.jld2", ind))
    end

    ϕ = load_object(filename)

    return ϕ
end

"""
    surfaces_from_file(dir::String)

Loads surfaces from file.
"""
function surfaces_from_file(dir::String)
    return load_object(dir)
end
