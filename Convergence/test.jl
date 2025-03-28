
#basic test on how finite difference error works.
using Plots; plotlyjs()
using Statistics
#%%


f(x) = sin(x)

df(x) = cos(x)


Df(x, h) = (f(x+h) - f(x-h)) / (2*h)
Df4(x, h) = (-f(x+2*h) + 8 * f(x+h) - 8 *f(x-h) + f(x-2*h)) / (12*h)

xvals = LinRange(0, 2π, 100)
#%%

hlist = @. 0.1 ^ (1:10)
hlist = [0.1, 0.05, 0.03, 0.01, 0.005, 0.003, 0.001]
Nh = length(hlist)
error = zeros(length(hlist), length(xvals))


for (j, h) in enumerate(hlist)
        
    for (i, x) in enumerate(xvals)

        error[j, i] = abs(Df4(x, h) - df(x))
    end
end
#%%
avg_error = zeros(Nh)
for i in 1:Nh
    avg_error[i] = mean(error[i, :])
end

#%%
lh = log.(hlist)
le = log.(avg_error)
plot(-log.(hlist), -log.(avg_error))
ind = 2
display((le[ind+1] - le[ind])/(lh[ind+1] - lh[ind]))
grad = zeros(Nh - 1)
for i in 1:Nh - 1
    grad[i] = (le[i+1] - le[i])/(lh[i+1] - lh[i])
end

plot(-lh[1:end-1], grad)
