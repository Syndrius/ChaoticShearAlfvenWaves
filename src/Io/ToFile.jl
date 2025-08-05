"""
    inputs_to_file(; prob::ProblemT, grids::GridsT, dir::String)

Writes the prob and grids structs to file using JLD2.
"""
function inputs_to_file(; prob::ProblemT, grids::GridsT, solver::SolverT, dir::String)

    save_object(dir*"prob.jld2", prob)
    save_object(dir*"grids.jld2", grids)
    save_object(dir*"solver.jld2", solver)

end


"""
    evals_to_file(evals::EvalsT, dir::String)

Writes the evals struct to file using JLD2.
"""
function evals_to_file(evals::EvalsT, dir::String)

    save_object(dir*"evals.jld2", evals)

end

"""
    efuncs_to_file(ϕ, ϕft, dir::String)

Writes the eigenfunctions and the fourier transformed eigenfunctions
to file. Each eigenfunction is written to an individual file.
"""
function efuncs_to_file(ϕ::Array{ComplexF64, 4}, ϕft::Array{ComplexF64, 4}, dir::String)
    nevals = size(ϕft)[1]
    efunc_base = dir*"/efuncs/"
    efunc_ft_base = dir*"/efuncs_ft/"
    mkpath(efunc_base)
    mkpath(efunc_ft_base)
    for i in 1:nevals

        efunc_label = efunc_base * @sprintf("efunc%05d.jld2", i)
        efunc_ft_label = efunc_ft_base * @sprintf("efunc%05d.jld2", i)

        if !isnothing(ϕ)
            #handles ϕ for fss case
            save_object(efunc_label, ϕ[i, :, :, :])
        end
        save_object(efunc_ft_label, ϕft[i, :, :, :])

    end
    
end

function efuncs_to_file(ϕ::Array{ComplexF64, 5}, ϕft::Array{ComplexF64, 5}, dir::String)
    nevals = size(ϕft)[1]
    efunc_base = dir*"/efuncs/"
    efunc_ft_base = dir*"/efuncs_ft/"
    mkpath(efunc_base)
    mkpath(efunc_ft_base)
    for i in 1:nevals

        efunc_label = efunc_base * @sprintf("efunc%05d.jld2", i)
        efunc_ft_label = efunc_ft_base * @sprintf("efunc%05d.jld2", i)

        if !isnothing(ϕ)
            #handles ϕ for fss case
            save_object(efunc_label, ϕ[i, :, :, :, :])
        end
        save_object(efunc_ft_label, ϕft[i, :, :, :, :])

    end
    
end
