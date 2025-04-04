
#fortran is succesfully able to solve and write data, but can we actually look at the damn stuff
using MID
using DelimitedFiles



dir_base = "/Users/matt/mid_data_test/"


fortran_process(dir_base)

efuncs = par_funcs_from_file("/Users/matt/code/efuncs.dat", 100, grids)


filename = "/Users/matt/code/efunc00003.txt"

nevals = 1
#actually using a petsc viewer now, currently with Ascii, which is terrible, will use hdf5 later

n = 8 * grids.r.N * grids.θ.N * grids.ζ.N

efuncs = open(filename, "r") do file
    nevals = 1
    efuncs = Array{ComplexF64}(undef, n, nevals) 
    #temp array for combining two floats into a single complex number.
    tmp = Array{ComplexF64}(undef, n) 
    s = readlines(file)
    #s = vcat(ts[1:723], ts[725:end])
    for i in 1:nevals

        j = 1
        
        #awkward pattern avoids Slepc header before each vector.
        #for (j, str) in enumerate(split.(s[4+(i-1)*(n+2):2*i+i*n]))
        for str in split.(s[4:end])

            if str[1] == "Process"
                display("found the process")
                #tmp[j] = 0.0 + 0.0im
                continue
            end
            
            #display(str)
            real = parse(Float64, str[1])
            imag = parse(Float64, str[3]) #* 1im
            #display(real)
            #display(imag)
            if str[2]=="+"
                tmp[j] = real + imag*1im
            else 
                tmp[j] = real - imag*1im
            end
            j += 1
        end
        efuncs[:, i] = tmp
    end

    return efuncs
end


ϕp = zeros(ComplexF64, grids.r.N, grids.θ.N, grids.ζ.N)
ϕpft = zeros(ComplexF64, grids.r.N, grids.θ.N, grids.ζ.N)

plan = MID.PostProcessing.create_ft_plan(ϕpft, grids)

#efunc = zeros(ComplexF64, 8 * 180)

#efunc[1:8:end] = efuncs[:, 1]

MID.PostProcessing.reconstruct_phi!(efuncs[:, 1], grids, ϕp, ϕpft, plan)


potential_plot(ϕpft, grids)



a = 1















ϕft=parse.(ComplexF64, open(readdlm,file)[1, :])

data = open(readdlm,file)


function pair_to_complex(pair)
    display(pair)
    #display((split(pair[1][2:end-1], ",")))
    nums = split(pair[1][2:end-1], ",")
    return parse(Float64, nums[1]) + parse(Float64, nums[2])*1im 
end
ϕ = map(pair_to_complex, (open(readdlm,file)[3:end], " "))

efunc = zeros(ComplexF64, 8 * 180)

efunc[1:8:end] = ϕ

ϕp = zeros(ComplexF64, grids.r.N, grids.θ.N, grids.ζ.N)
ϕpft = zeros(ComplexF64, grids.r.N, grids.θ.N, grids.ζ.N)

plan = MID.PostProcessing.create_ft_plan(ϕpft, grids)

MID.PostProcessing.reconstruct_phi!(ϕ, grids, ϕp, ϕpft, plan)


potential_plot(ϕpft, grids)
