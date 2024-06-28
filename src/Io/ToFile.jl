
"""
    problem_to_file(; prob::ProblemT, filename::String)

Write the problem struct to file.
"""
function problem_to_file(; prob::ProblemT, filename::String)

    open(filename, "w") do file
        write(file, "Problem:\n")
        write(file, "q profile: " * string(prob.q) * "\n")
        write(file, "metric: " * string(prob.compute_met) * "\n")
        write(file, "density: " * string(prob.dens) * "\n")
        write(file, "Island:\n")
        write(file, " - m0: " * string(prob.isl.m0) * "\n")
        write(file, " - n0: " * string(prob.isl.n0) * "\n")
        write(file, " - A: " * string(prob.isl.A) * "\n")
        write(file, "Geometry:\n")
        write(file, " - R0: " * string(prob.geo.R0) * "\n")
        write(file, " - a: " * string(prob.geo.a) * "\n")
        write(file, " - B0: " * string(prob.geo.B0) * "\n")
        write(file, "Resisitivity: " * string(prob.δ))
    end

end


"""
    grids_to_file(; grids::GridsT, filename::String)

Write the grids struct to file.
"""
function grids_to_file(; grids::GridsT, filename::String)
    #maybe we should force this to be txt file?
    open(filename, "w") do file
        if grids.θ isa SMGridDataT
            write(file, "Grids: FSS\n")
        else
            write(file, "Grids: FFS\n")
        end
        write(file, "Radial Grid:\n")
        grid_to_file(file, grids.r)
        write(file, "Poloidal Grid:\n")
        grid_to_file(file, grids.θ)
        write(file, "Toroidal Grid:\n")
        grid_to_file(file, grids.ζ)
        
    end
end

"""
    grid_to_file(file, grid::FEMGridDataT)

Writes a finite element method grid to file.
"""
function grid_to_file(file, grid::FEMGridDataT)

    write(file, " - N: " * string(grid.N) * "\n")
    write(file, " - sep1: " * string(grid.sep1) * "\n")
    write(file, " - sep2: " * string(grid.sep2) * "\n")
    write(file, " - frac: " * string(grid.frac) * "\n")
    write(file, " - gp: " * string(grid.gp) * "\n")
    write(file, " - pf: " * string(grid.pf) * "\n")

end

"""
    grid_to_file(file, grid::SMGridDataT)

Writes a spectral method grid to file.
"""
function grid_to_file(file, grid::SMGridDataT)

    write(file, " - start: " * string(grid.start) * "\n")
    write(file, " - count: " * string(grid.count) * "\n")
    write(file, " - incr: " * string(grid.incr) * "\n")
    write(file, " - f_quad: " * string(grid.f_quad) * "\n")

end



"""
    inputs_to_file(; prob::ProblemT, grids::GridsT, dir::String

Writes the problem and grids structs to files in the given directory.
"""
function inputs_to_file(; prob::ProblemT, grids::GridsT, dir::String)

    problem_to_file(prob=prob, filename=dir*"prob.txt")
    grids_to_file(grids=grids, filename=dir*"grids.txt")
end


"""
    eigvals_to_file(; ω::Array{ComplexF64}, filename::String)

Writes the eigenvalues to file.
"""
function eigvals_to_file(; ω::Array{ComplexF64}, filename::String)
    open(filename, "w") do file
        writedlm(file, ω, ",")
    end
end

"""
    eigfuncs_to_file(; ϕ::Array{ComplexF64, 2}, filename::String)

Writes the eigenfunctions to file.
"""
function eigfuncs_to_file(; ϕ::Array{ComplexF64, 2}, filename::String)
    open(filename, "w") do file
        writedlm(file, ϕ, ",")
    end
end