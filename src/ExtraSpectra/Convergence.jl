"""
Runs a convergence test, where different N and δ values are considered and results are written to file.
This assumes only N and δ are varied.

# Args
grids::GridsT - Grids used, will be modified with each N value.
prob::ProblemT - Problem to solve, will be modified with each δ value.
Nlist::Array{Int64} - List of N's to test with.
δlist::Array{Float64} - List of δ's to test with.
dir_base::String - Directory to save data to.
σ::Float64 - Target frequency.
"""
function convergence_test(; grids::GridsT, prob::ProblemT, Nlist::Array{Int64}, δlist::Array{Float64}, dir_base::String, σ::Float64)

    #assumes only N and δ will change

    mkpath(dir_base) #should make the directory and not poo itself if it already exists.

    for N in Nlist
        
        #@reset allows immutable structs to be redefined.
        @reset grids.rd.N = N #can we do this?
        for δ in δlist

            @reset prob.δ = δ

            #extracts the power of δ to use as a label.
            dlab = parse(Int64, split(string(δ), "-")[end])
    
            dir = dir_base * @sprintf("N_%s_delta_%s/", N, dlab)
            mkdir(dir) #maybe this should be mkpath, but mkdir seems better in this case.
    
            inputs_to_file(prob=prob, grids=grids, dir=dir)
    
            ω, ϕ = construct_and_solve(prob=prob, grids=grids, full_spectrum=false, σ=σ, reconstruct=false);
    
            eigvals_to_file(ω=ω, filename=dir*"vals.dat")
    
            eigfuncs_to_file(ϕ=ϕ, filename=dir*"funcs.dat")
    
        end
    end

end


"""
Reads the output data from the convergence_test function. Assumes only the first eigenvalue is relevant.

# Args
Nlist::Array{Int64} - List of N's used.
δlist::Array{Float64} - List of δ's used.
dir_base::String - Directory where data was saved.
"""
function read_convergence_data(; Nlist::Array{Int64}, δlist::Array{Float64}, dir_base::String)

    ωlist = zeros(ComplexF64, length(Nlist), length(δlist))

    for (j, N) in enumerate(Nlist)

        for (i, δ) in enumerate(δlist)

            dlab = parse(Int64, split(string(δ), "-")[end])
            dir = dir_base * @sprintf("N_%s_delta_%s/", N, dlab)
    
            vals_file = dir * "vals.dat"
            #funcs_file = dir * "funcs.dat"
            
            ω = readdlm(vals_file, ',', ComplexF64)
            #doesn't do this, ideally we probably want some function that can read the ith eigenfunction only.
            #ϕ = reconstruct_phi(readdlm(funcs_file, ',', ComplexF64), length(ω), N, grids.pmd.count, grids.tmd.count)
    
            #most of these are the first index, some weird stuff is going on later though!
            ωlist[j, i] = ω[1]
            
            #doesn't do this here!
            #plot_potential(r=rgrid, ϕ=ϕ, ind=1, pmd=grids.pmd, n=1)
    
        end
    end

    return ωlist
end


"""
Same as convergence_test() but specifically for the two mode case, used for comparison. Defaults to Singular Bowden case.

# Args
Nlist::Array{Int64} - List of N's to test with.
δlist::Array{Float64} - List of δ's to test with.
dir_base::String - Directory to save data to.
R0::Float64 - Major radius.
σ::Float64 - Target frequency.
sep1::Float64=0.75 - Start of clustered region. 
sep2::Float64=0.8 - End of clustered region.
frac::Float64=0.2 - Proportion of clustered region. 
m::Int64=1 - First poloidal mode.
n::Int64=-1 - Toroidal mode.
"""
function two_mode_convergence(; Nlist::Array{Int64}, δlist::Array{Float64}, dir_base::String, R0::Float64, σ::Float64, sep1::Float64=0.75, sep2::Float64=0.8, frac::Float64=0.2, m::Int64=1, n::Int64=-1)

    mkpath(dir_base)

    for N in Nlist

        rgrid = clustered_grid(N, 0.75, 0.8, 0.2);

        for δ in δlist

            dlab = parse(Int64, split(string(δ), "-")[end])
            dir = dir_base * @sprintf("N_%s_delta_%s/", N, dlab)

            mkdir(dir) 

            ω, _, _ = two_mode(rgrid=rgrid, R0=R0, m=m, n=n, δ=δ, σ=σ);

            eigvals_to_file(ω=ω, filename=dir*"vals.dat")
        end
    end
end
