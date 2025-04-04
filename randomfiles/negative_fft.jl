
#wot the fek is going on with the fft.
using FFTW

N = 10
#so this N+1 [1:end-1] is v essential for clean ft. May need to double check in our case.
x = LinRange(0, 2π, N+1)[1:end-1]
y = LinRange(0, 2π, N+1)[1:end-1]

f(x, y) = (cos(4*x) + cos(3*x) + cos(x)*cos(y) - cos(2*y) + cos(3*y)) * exp(1im*x)#+1 #+ cos(2*x)*cos(2*y*2x)

g = zeros(ComplexF64, N, N)
for i in 1:N, j in 1:N
    g[i, j] = f(x[i], y[j])
end


contourf(x, y, g')


res = fft(g, [1, 2]);


show(stdout, "text/plain", real.(round.(res)))
#contourf(x, y, real(res))

scatter(0:N-1, real.(res[:, 2]), legend=false)

scatter(0:N-1, real.(res[13, :]), legend=false)




f(x) = exp(1im * x) + exp(1im * 4*x)
#cos = 1/2(exp + exp) so we get the negative signal as well, at half strength compared to exp.
f(x) = cos(x) + cos(4*x)


res = fft(f.(x))


#display(res)

scatter(0:N-1, real.(res))
plot!(x, imag.(res))
println(res)


#NOW inverse fft

#say we have 
mlist = [1, 2, 3, -1, -2]
#with 
vals = [1.0, 0.0, 0.0, 0.0, 0.0]
#need to create array of [0, 1, 2, 3, -3, -2, -1]
maxm = maximum(abs.(mlist))

maxm = 100

grid = LinRange(0, 2π, 2*maxm+1+1)[1:end-1]


res = zeros(2*maxm+1)

for i in 1:1:length(mlist)

    res += vals[i] * exp.(1im * mlist[i] .* grid)
end

xplot = LinRange(0, 2π, length(res)+1)[1:end-1]
plot(xplot, real.(res))
plot(xplot, imag.(res))





fvals = zeros(2*maxm+1)

for i in 1:1:length(mlist)

    ind = maxm + mlist[i] + 1
    fvals[ind] = vals[i]
end
println(fvals)


res = ifft(ifftshift(fvals))
xplot = LinRange(0, 2π, length(res)+1)[1:end-1]
plot(xplot, real.(res))
plot(xplot, imag.(res))


maxm = 100

#won't work with incr≠1, can probs use incr in here though???
#may cause issues with zeros? I guess we can just make a huge array and set heaps of values to zero.
fms = collect(-maxm:maxm)
fms = vcat(0:maxm, -maxm:-1)
fvals = zeros(length(fms))

for i in 1:1:length(mlist)
    find = mode_to_array(mlist[i], fms)
    display(find)
    fvals[find] = vals[i]

end
println(fvals)

#not an ideal way of doing this but I think it will work.
function mode_to_array(mode, fms)
    #display(mode)
    if mode >= 0
        return fms[mode+1] + 1
    else
        display(fms[end+mode+1])
        return fms[end+mode+1] + 1
    end
end

res = ifft(fvals)
xplot = LinRange(0, 2π, length(res)+1)[1:end-1]
plot(xplot, real.(res))


#say we have array 
N = 100
x = zeros(N)

x[3] = N
x[end-1] = N

y = ifft(x)

xplot = LinRange(0, 2π, N+1)[1:end-1]
plot(xplot, real.(y))
#plot(xplot, imag.(y))
