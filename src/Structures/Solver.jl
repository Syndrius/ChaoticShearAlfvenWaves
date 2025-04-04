"""
Struct storing the solving parameters for shift and invert.

### Fields
- nev::Int64=100 - Number of eigenvalues to find.
- target::Float64=0.0 - Target frequency. Solution will be the nev nearest eigenvalues to this target.
"""
@kwdef struct ShiftInvertSolverT <: SolverT
    nev :: Int64 = 100
    target :: Float64 = 0.0
end


"""
Struct storing the solving parameters for slice solving. 

### Fields
- nev::Int64=100 - Number of eigenvalues to find for each slice.
- targets::Array{Float64} - List of the target frequency. Solution will be the nev nearest eigenvalues to each target in targets.
"""
@kwdef struct SliceSolverT <: SolverT
    nev :: Int64 = 100
    targets :: Array{Float64}
end


"""
Struct storing the solving parameters.

### Fields
- ideal::Bool=false - Boolean for the problem being ideal or not meaning the matrices are Hermitian or not.
"""
@kwdef struct FullSpectrumSolverT <: SolverT
    ideal :: Bool = false
end


"""
    init_solver(;nev :: Int64 = 100, target :: Float64 = 0.0, nshifts::Int64=1, tlow::Float64=0.0, thigh::Float64=1.0, targets::Array{Float64}=Float64[], full_spectrum::Bool=false, prob::ProblemT)

Constructor for struct which stores the key information that defines how the problem will be solved.

# Args
- nev::Int64=100 - Number of eigenvalues to solve.
- target::Float64=0.0 - target for the shift and invert transformation.
- nshifts::Int64=1 - Number of slices to take for slice solving. If != 1, a list of targets from between tlow and thigh is created.
- tlow::Float64=0.0 - Lowest target for slice list.
- thigh::Float64=1.0 - Highest target for slice list.
- targets::Array{Float64}=Float64[] - List of targets of slice solving.
- full_spectrum::Bool=false - Set true to solve the full spectrum, only practical for small grids.
- prob::ProblemT - problem struct for normalising and determining if matrices are Hermitian.
"""
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


"""
    get_targets(nshifts::Int64, tlow::Float64, thigh::Float64, geo::GeoParamsT)

Function that creates the target list for slice solving.
"""
function get_targets(nshifts::Int64, tlow::Float64, thigh::Float64, geo::GeoParamsT)
    #think just passing an array of targets is almost always the better option tbh
    #but good to have the option.

    #not certain that we do want to do this tbh.
    #tlow = tlow^2 / geo.R0^2
    #thigh = thigh^2 / geo.R0^2

    return LinRange(tlow, thigh, nshifts) .^2 ./ geo.R0^2
end
