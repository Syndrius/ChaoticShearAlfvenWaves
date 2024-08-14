
function inputs_from_file(; dir::String)

    grids = load_object(dir*"grids.jld2")
    prob = load_object(dir*"prob.jld2")
    
    return prob, grids

end

function evals_from_file(; dir::String)

    evals = load_object(dir*"evals.jld2")

    return evals
end


function efunc_from_file(; dir::String, ind, ft=true)

    if ft
        filename = dir * @sprintf("efuncs_ft/efunc%04d.jld2", ind)
    else
        filename = dir * @sprintf("efuncs/efunc%04d.jld2", ind)
    end

    ϕ = load_object(filename)

    return ϕ
end