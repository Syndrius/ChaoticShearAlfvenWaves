#understanding fft integration, making sure we have the correct constants scaling our terms in the qfm calculations
using FFTW
using Plots; plotlyjs()
#%%
N = 100
a = 3
x = LinRange(0, 2π*a, N+1)[1:end-1]

y = x .^2 .+ 1;

fy = rfft(y) ./ N .* 2*π*a;
#%%

plot(0:(N ÷2), real.(fy))
plot(0:N÷2, imag.(fy))

#%%
#why is fft so fkn confusing

amp = @. (1-x)^2

y = @. amp' * exp(1im*2 * x);

fy = fft(y, [1]) ./ N

plot(x, real(y))
#%%
p = plot()
for i in 1:N
    plot!(x, real(fy[i, :]))
end
display(p)

