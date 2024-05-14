


#reading and writing problems and grids to files.
#we have chosen to do this quite verbosly and manually, making the txt files much easier to read. But any slight changes will completly cook everything!

function problem_to_file(; prob::ProblemT, filename::String)

    #qname = 
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

function grids_to_file(; grids::GridsT, filename::String)
    open(filename, "w") do file
        write(file, "Grids:\n")
        write(file, "Radial Grid:\n")
        #this may get stpuid as the grid gets big, oh well.
        write(file, " - grid: " * string(grids.rd.grid) * "\n")
        write(file, " - gp: " * string(grids.rd.gp) * "\n")
        write(file, " - N: " * string(grids.rd.N) * "\n")
        write(file, "Poloidal Grid:\n")
        write(file, " - start: " * string(grids.pmd.start) * "\n")
        write(file, " - count: " * string(grids.pmd.count) * "\n")
        write(file, " - incr: " * string(grids.pmd.incr) * "\n")
        write(file, " - f_quad: " * string(grids.pmd.f_quad) * "\n")
        write(file, "Toroidal Grid:\n")
        write(file, " - start: " * string(grids.tmd.start) * "\n")
        write(file, " - count: " * string(grids.tmd.count) * "\n")
        write(file, " - incr: " * string(grids.tmd.incr) * "\n")
        write(file, " - f_quad: " * string(grids.tmd.f_quad) * "\n")
    end
end



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

function grids_from_file(; filename::String)
    file = open(filename, "r") 
    readline(file) #ignore first line

    #Radial Grid
    readline(file)
    s = readline(file)
    #display(s)
    rgrid = parse.(Float64, split(chop(s[10:end]; head=1, tail=1), ','))
    #rgrid = parse(Array{Float64}, s[10:end])
    s = readline(file)
    #display(s)
    gp = parse(Int64, s[8:end])
    s = readline(file)
    N = parse(Int64, s[7:end])

    #Poloidal Grid
    readline(file)
    s = readline(file)
    pstart = parse(Int64, s[11:end])
    s = readline(file)
    pcount = parse(Int64, s[11:end])
    s = readline(file)
    pincr = parse(Int64, s[10:end])
    s = readline(file)
    pf_quad = parse(Int64, s[12:end])

    #Toroidal Grid
    readline(file)
    s = readline(file)
    tstart = parse(Int64, s[11:end])
    s = readline(file)
    tcount = parse(Int64, s[11:end])
    s = readline(file)
    tincr = parse(Int64, s[10:end])
    s = readline(file)
    tf_quad = parse(Int64, s[12:end])

    rd = RDataT(grid=rgrid, gp=gp, N=N)
    pmd = ModeDataT(start=pstart, count=pcount, incr=pincr, f_quad=pf_quad)
    tmd = ModeDataT(start=tstart, count=tcount, incr=tincr, f_quad=tf_quad)

    return GridsT(rd = rd, pmd=pmd, tmd=tmd)
end


#these seem to be working!
function inputs_from_file(; dir::String)
    prob = problem_from_file(filename=dir * "prob.txt")
    grids = grids_from_file(filename=dir * "grids.txt")

    return prob, grids

end

function inputs_to_file(; prob::ProblemT, grids::GridsT, dir::String)

    #I guess we are assuming that dir ends with /
    problem_to_file(prob=prob, filename=dir*"prob.txt")
    grids_to_file(grids=grids, filename=dir*"grids.txt")
end


#straightforward write to file function
function eigvals_to_file(; ω, filename)
    open(filename, "w") do file
        writedlm(file, ω, ",")
    end
end

#straightforward write to file function
#may want to split this into real and imag to match parallel case.
function eigfuncs_to_file(; ϕ, filename)
    open(filename, "w") do file
        writedlm(file, ϕ, ",")
    end
end