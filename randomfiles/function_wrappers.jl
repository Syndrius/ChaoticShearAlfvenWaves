
#testing some things with function wrapper.s

using FunctionWrappers
using JLD2

#%%
struct test
    q :: FunctionWrappers.FunctionWrapper{Tuple{Float64, Float64}, Tuple{Float64}}
end

#%%
function test_q(r::Float64)

    return 1.0 + r^2, 2*r
end
function more_q(r::Float64, b::Float64)
    return b + r^2, 2*r
end
#%%
mystruct = test(test_q)

anon_q(r) = more_q(r, 2.0)

my_s = test(anon_q)

@code_warntype mystruct.q(0.5)

mystruct.q(0.5)
my_s.q(0.5)

save_object("test_struct.jld2", mystruct)
save_object("testy_struct.jld2", my_s)

#this perhaps exagerates our problemo
#mayhap a blessing in disguise?
l_struct = load_object("test_struct.jld2")
l_s = load_object("testy_struct.jld2")

l_struct.q(0.5)
l_s.q(0.5)

#%%
function test_met!(s::Float64)
    display(2*s)
end

struct met_struct3
    mm :: FunctionWrappers.FunctionWrapper{Nothing, Tuple{Float64}}
end
#%%

st = met_struct3(test_met!)

st.mm(1.6)
