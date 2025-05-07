#trying to understand real fft better.
using FFTW
using Plots; plotlyjs()

#%%
x = LinRange(0, 2π, 100)

y1 = @. sin(x) + 3*sin(4*x) + exp(im * x) + 0.01;
y = @. sin(2*x) + 3*sin(-5*x) + 1*sin(1*x) + 0.01;
y1 = @. exp(im * 2*x) + 3*exp(-im*5*x) + exp(im * x) + 0.01;

#%%

fty = rfft(y);
fty1 = fft(y);

length(fty)
length(fty1)

plot(x, real.(y))
plot(x, imag.(y))
plot(0:99, real.(fty1), legend=false)
plot(0:50, real.(fty)) 
plot(0:99, imag.(fty1))


res = irfft(fty, 100);
res1 = irfft(fty, 101);

plot(res)
plot(res1)

length(fty)
length(res)
