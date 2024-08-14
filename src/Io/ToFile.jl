
function inputs_to_file(; prob::ProblemT, grids::GridsT, dir::String)

    save_object(dir*"prob.jld2", prob)
    save_object(dir*"grids.jld2", grids)

end

function evals_to_file(evals::EvalsT, dir::String)

    save_object(dir*"evals.jld2", evals)

end


function efuncs_to_file(ϕ, ϕft, dir::String)
    nevals = size(ϕft)[1]
    efunc_base = dir*"/efuncs/"
    efunc_ft_base = dir*"/efuncs_ft/"
    mkpath(efunc_base)
    mkpath(efunc_ft_base)
    for i in 1:nevals

        efunc_label = efunc_base * @sprintf("efunc%04d.jld2", i)
        efunc_ft_label = efunc_ft_base * @sprintf("efunc%04d.jld2", i)

        if !isnothing(ϕ)
            #handles ϕ for fss case
            save_object(efunc_label, ϕ[i, :, :, :])
        end
        save_object(efunc_ft_label, ϕft[i, :, :, :])

    end
    
end