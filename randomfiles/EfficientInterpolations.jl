#our interpolation is slow af.
#here we try and make it faster, by creating the full interpolation object first
#Julia packages tend to create structs and use something like
#function (a::structT)(x) 

#don't think this will actually help us though!
#think this won't be worth the hassle

#again we will start in 1d
N = 25
xvals = LinRange(0, 2π, N+1)[1:end-1]
f(x) = sin(x)
df(x) = cos(x)
fvals = zeros(N, 2)
fvals[:, 1] = f.(xvals)
fvals[:, 2] = df.(xvals)
#%%


#%%

function h00(t)

    return 2*t^3 - 3*t^2 + 1
    #return (1+2*t)*(1-t)^2
end

function h10(t)
    #return 2 * (t^3-2*t^2+t)
    return t*(1-t)^2
end

function h01(t)
    return -2t^3+3t^2
    #return t^2*(3-2*t)
end

function h11(t)
    #return 2*(t^3-t^2)
    return t^2*(t-1)
end

function hb(t, h, dt)
    #additional jacobian term is very awkward.

    if h==1
        return h00(t)
    elseif h==2
        return h10(t) * dt
    elseif h==3
        return h01(t)
    else
        return h11(t) * dt
    end
end
#%%
function create_hi(fvals, xvals)

    #this seems mad as well.
    hi = Array{Function}(undef, length(xvals))

    for i in 1:length(xvals)-1
        dx = xvals[i+1] - xvals[i]
        #f(x) = fvals[i, 1] * h00(x) + fvals[i, 2] * h10(x) * dx + fvals[i+1, 1] * h01(x) + fvals[i+1, 2] * h11(x) * dx
        (hi[i])(x) = fvals[i, 1] * h00(x) + fvals[i, 2] * h10(x) * dx + fvals[i+1, 1] * h01(x) + fvals[i+1, 2] * h11(x) * dx
        #hi[i] = f
    end
    dx = 2π + (xvals[1] - xvals[end])
    f(x) = fvals[end, 1] * h00((x-xvals[end])/dx) + fvals[end, 2] * h10((x-xvals[end])/dx) * dx + fvals[1, 1] * h01((x-xvals[end])/dx) + fvals[1, 2] * h11((x-xvals[end])/dx) * dx
    f(x) = 2*x^10
    hi[end] = f
    
    return hi
end
        

hi = create_hi(fvals, xvals)
#%%

function eval_hi(hi, x, xvals)

    ind = argmin(abs.(x .- xvals))
    display(ind)
    if x >= xvals[ind]
        node = ind
    else
        node = ind-1
    end
    t = (x - xvals[node]) / (xvals[node+1] - xvals[node-1])
    return hi[node](t)
end
eval_hi(hi, 0.3, xvals)
display(xvals)

f(0.3)
#%%
Ni = 200
xi = LinRange(0.5, π, Ni+1)[1:end-1]
ϕ = zeros(Ni)
fi = f.(xi)
for (i, x) in enumerate(xi)
    ϕ[i] = eval_hi(hi, x, xvals)
end
#%%
plot(xi, fi)
plot!(xi, ϕ)
