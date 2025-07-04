#trying to understand real fft better.
using FFTW
using Plots; plotlyjs()

#%%
N = 50
a = 3
x = LinRange(0, 2*a*π, N)

y1 = @. sin(x/a) + 3*sin(4*x/a) + exp(im * x/a) + 0.01;
y1 = @. exp(im * 2*x/a) + 3*exp(-im*5*x/a) + exp(im * x/a) + 0.01;

y = @. 10*sin(2*x/a) + 3*sin(-5*x/a) + 1*sin(1*x/a) + 0.01; #+ 1im*cos(2*x/a);
#%%
fty1 = fft(y);
fty = rfft(y); #so if the input is purely real, the spectrum is symmetric, so we can ignore half.
plot(0:N-1, real.(fty1), legend=false)
plot(0:floor(N/2), real.(fty), legend=false)
#%%
fake_fft = [fty[1:end-1]; fty[end:-1:2]];
length(fake_fft)
plot(0:N-1, real.(fake_fft))
plot(0:N-1, real.(fty1))
display(real.(fake_fft) .- real.(fty1))
#%%
#now the inverse
yrinv = irfft(fty, N);
yinv = ifft(fty1);

length(yinv) #N
length(yrinv) #N

plot(x, yrinv)
plot(x, real.(yinv))

plot(x, y)


#%%
fty1 = fft(y);

length(fty) #N/2 + 1
length(fty1) #N

plot(x, real.(y))
plot(x, imag.(y))
plot(0:N-1, real.(fty1), legend=false)
plot(0:50, real.(fty)) 
plot(0:99, imag.(fty1))


res = irfft(fty, 100);
res1 = irfft(fty, 101);

plot(res)
plot(res1)

length(fty)
length(res)
