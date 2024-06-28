


#reading and writing problems and grids to files.
#we have chosen to do this quite verbosly and manually, making the txt files much easier to read. But any slight changes will completly cook everything!



"""
    problem_from_file(; filename::String)

Reads the problem struct from file.
"""
function problem_from_file(; filename::String)

    file = open(filename, "r") 
    readline(file) #ignore first line
    s = readline(file)
    q_prof = getfield(Main, Symbol(s[12:end])) #not ideal!
    s = readline(file)
    met = getfield(Main, Symbol(s[9:end]))

    s = readline(file)
    dens = getfield(Main, Symbol(s[10:end]))


    #island
    readline(file) #ignore
    s = readline(file)
    m0 = parse(Int64, s[8:end])
    s = readline(file)
    n0 = parse(Int64, s[8:end])
    s = readline(file)
    A = parse(Float64, s[7:end])

    #geometry
    readline(file)
    s = readline(file)
    R0 = parse(Float64, s[8:end])
    s = readline(file)
    a = parse(Float64, s[7:end])
    s = readline(file)
    B0 = parse(Float64, s[8:end])

    s = readline(file)
    δ = parse(Float64, s[15:end])
    close(file)

    isl = IslandT(A=A, m0=m0, n0=n0)
    geo = GeoParamsT(R0=R0, a=a, B0=B0) #a and B0 will not do anything!

    return ProblemT(q=q_prof, compute_met=met, dens=dens, isl=isl, geo=geo, δ=δ)
            
end


"""
    grids_from_file(; filename::String)

Reads the grids struct from file.
"""
function grids_from_file(; filename::String)
    file = open(filename, "r") 
    s = readline(file) 
    type = split(split(s, " ")[2], "")

    #Radial Grid:
    readline(file) 
    rgrid = grid_from_file(file, type[1])

    #Poloidal Grid
    readline(file) 
    θgrid = grid_from_file(file, type[2])
    
    #Toroidal Grid
    readline(file)
    ζgrid = grid_from_file(file, type[3])
    

    return init_grids(rgrid, θgrid, ζgrid)
end

"""
    grid_from_file(file::IOStream, type::SubString)

Reads a grid from file, based on type = "F" or "S".
"""
function grid_from_file(file::IOStream, type::SubString)

    if type == "F" || type == 'F'

        s = readline(file)
        N = parse(Int64, s[7:end])
        s = readline(file)
        sep1 = parse(Float64, s[10:end])
        s = readline(file)
        sep2 = parse(Float64, s[10:end])
        s = readline(file)
        frac = parse(Float64, s[10:end])
        s = readline(file)
        gp = parse(Int64, s[8:end])
        s = readline(file)
        pf = parse(Int64, s[8:end])

        return FEMGridDataT(N, sep1, sep2, frac, gp, pf)
    else
        s = readline(file)
        start = parse(Int64, s[11:end])
        s = readline(file)
        count = parse(Int64, s[11:end])
        s = readline(file)
        incr = parse(Int64, s[10:end])
        s = readline(file)
        f_quad = parse(Int64, s[12:end]) 
        return SMGridDataT(start, count, incr, f_quad)
    end
end



"""
    inputs_from_file(; dir::String)

Reads the problem and grids structs from files in the directory.
"""
function inputs_from_file(; dir::String)
    prob = problem_from_file(filename=dir * "prob.txt")
    grids = grids_from_file(filename=dir * "grids.txt")
    return prob, grids
end