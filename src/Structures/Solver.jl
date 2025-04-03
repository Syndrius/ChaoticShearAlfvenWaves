#think we should get this stuff working for shiftinvert solve in serial first.
#may help make the parallel version work better.
#so I think we will offer 3 (4) ways to solve. 
#in serial we also offer the ability to solve the full spectrum
#solve for a single target
#solve over a range with nshift targets inbetween tlow and thigh
#solve for a list of targets
#this will need work unfortuantly.

@kwdef struct ShiftInvertSolverT <: SolverT
    nev :: Int64 = 100
    target :: Float64 = 0.0
end

@kwdef struct SliceSolverT <: SolverT
    nev :: Int64 = 100
    targets :: Array{Float64}
end

@kwdef struct FullSpectrumSolverT <: SolverT
    ideal :: Bool = false
end


function init_solver(;nev :: Int64 = 100, target :: Float64 = 0.0, nshifts::Int64=1, tlow::Float64=0.0, thigh::Float64=1.0, targets::Array{Float64}=Float64[], full_spectrum::Bool=false, prob::ProblemT)

    if full_spectrum
        if prob.flr.δ == 0.0 && prob.flr.ρ_i == 0 && prob.flr.δ_e == 0
            return FullSpectrumSolverT(ideal=true)
        else
            return FullSpectrumSolverT(ideal=false)
        end
    elseif isempty(targets) && nshifts==1
        target = target^2 / prob.geo.R0^2
        return ShiftInvertSolverT(nev, target)
    elseif !isempty(targets)
        targets = @. targets^2 / prob.geo.R0^2
        return SliceSolverT(nev, targets)
    else
        targets = get_targets(nshifts, tlow, thigh, prob.geo)
        return SliceSolverT(nev, targets)
    end
end


function get_targets(nshifts::Int64, tlow::Float64, thigh::Float64, geo::GeoParamsT)
    #think just passing an array of targets is almost always the better option tbh
    #but good to have the option.

    #not certain that we do want to do this tbh.
    #tlow = tlow^2 / geo.R0^2
    #thigh = thigh^2 / geo.R0^2

    return LinRange(tlow, thigh, nshifts) .^2 ./ geo.R0^2
end
