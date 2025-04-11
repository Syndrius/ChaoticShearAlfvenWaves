#file for storing the collection of real fft function we use.
#bit of a disaster, doesn't help that many of the arrays are an unknown size
#and some of the ft are a complete mystery.
#probably gotta be a problemo for another day!


#takes the real fft in 1D.
function rfft1D!(cos_coef::Array{Float64}, sin_coef::Array{Float64}, func::Array{Float64}, fft_res::Array{ComplexF64}, plan::FFTW.rFFTWPlan)

    #shouldn't be doing this.
    Nfft = length(func)


    fft_res .= plan * func

    cos_coef .= real.(fft_res) ./ Nfft .* 2

    cos_coef[1] /= 2

    sin_coef .= @. -imag(fft_res) / Nfft * 2

    sin_coef[1] = 0.0

end


#obvs this should not exist, it is still being used by grad_jm!
function rfft1D_simple(f)

    Nfft = size(f)[end]
    res = rfft(f)

    cosout = real.(res) ./ Nfft .* 2
    sinout = -imag.(res) ./ Nfft .* 2

    cosout[1] /= 2
    sinout[1] = 0.0

    return cosout, sinout
end

#reconstructs r and t from the combined fourier coefficients.
#takes in inverse real fourier transform to get r, θ from the coefficients.
#somewhat inconsistently named tbh. May be worth changing this or the other stuff.
function get_r_t!(r::Array{Float64}, θ::Array{Float64}, coefs::CoefficientsT, ift_plan::AbstractFFTs.ScaledPlan, ft_r1D::Array{ComplexF64}, ft_θ1D::Array{ComplexF64}, Ntor::Int64)

    #2 is an optional arg.
    Nfft = 2 * (Ntor - 1)

    #note that these arrays have been pre-padded with zeros to make sure the ifft works.
    ft_r1D[1:Ntor] .= @. (coefs.rcos - 1im * coefs.rsin) * Nfft
    ft_θ1D[1:Ntor] .= @. (coefs.θcos - 1im * coefs.θsin) * Nfft

    #unclear!
    ft_r1D[1] *= 2
    ft_θ1D[1] *= 2

    #take the inverse ft.
    r .= ift_plan * ft_r1D
    θ .= ift_plan * ft_θ1D
    
end

#still being used!
function rfft1D_JM(f)


    Nfft = size(f)[end]
    ffft = rfft(f, [2]) #for the JM case we need to do the second dim.
    #unsure why or if we can combine these functions.

    #println(ffft)

    cos_out = real.(ffft) ./ Nfft .* 2

    
    #1d case
    cos_out[:, 1] ./= 2
    

    sin_out = @. -imag(ffft) / Nfft * 2

    sin_out[:, 1] .= 0.0

    #print(sin_out)
    #probably a bad idea to have this here, but it is needed.
    return transpose(cos_out), transpose(sin_out)


end

function irfft1D(cos_in::Array{Float64, 1}, sin_in::Array{Float64, 1}, nfft_multiplier=2)


    Nfft = nfft_multiplier * (size(cos_in)[end] - 1)
    #display(Nfft * 2)
    ffft = (cos_in .- 1im .* sin_in) .* Nfft

    
    #display(Nfft)
    #display(length(ffft))

    
    #may need to change this
    #this is differretn, but I think because this is 1d it is ok.
    #will be a problemo later.
    #also do not understand why this is occuring.
    ffft[1] *= 2

    #println(ffft)
    #probbably need to pick the dim properly. Not sure what shape cos_in and sin_in have

    #need to pad ffft with zeros as it is not long enough for the inverse fft of the desired length
    #python does this automatically!
    Npad = floor(Int64,  Nfft / 2 - length(ffft)/2) + 1
    #ft_pad = [zeros(Npad); ffft; zeros(Npad)];

    #ok so this matches python, bit odd tbh.
    #so 2*Npad ; ffft with ifftshift makes this match python
    #I am sure there are other combos that would work.
    #so does ffft; 2*Npad, without anyshift.
    #pretty hekin annoying. But atleast we are there now!
    ft_pad = [ffft; zeros(Npad); zeros(Npad)];
    #may need to shift this!
    #display(length(ft_pad))
    #println(irfft(fftshift(ft_pad), 2 * Nfft) .* (Nfft / length(ffft)))
    #println(irfft(ifftshift(ft_pad), 2 * Nfft))# .* (Nfft / length(ffft)))
    #println(irfft(ft_pad, 2 * Nfft))# .* (Nfft / length
    #perhaps irfft is different in other ways???
    #display(length(ft_pad))
    #display(2*Nfft)
    return irfft(ft_pad, 2 * Nfft)

end

function rfft2D(f, mpol, ntor)

    #zhisongs code considers case where mpol and ntor are unspecified.

    Nfft1, Nfft2 = size(f)

    #going to have to assume a 2d function here.
    fftout = rfft(f, [1, 2])

    #display(fftout)

    fftcos = @. real(fftout) / Nfft1 / Nfft2 * 2

    fftcos[1, :] ./= 2
    fftcos[end, :] ./= 2

    fftsin = @. -imag(fftout) / Nfft1 / Nfft2 * 2
    fftsin[1, :] ./= 2
    fftsin[end, :] ./= 2

    #display(fftcos)

    mpol_data = div(Nfft1, 2)
    ntor_data = div(Nfft2, 2) - 1

    cn = zeros((mpol + 1, 2 * ntor + 1))
    sn = zeros((mpol + 1, 2 * ntor + 1))

    #oh so now we aren't doing any padding... cooked.
    if mpol < mpol_data && ntor < ntor_data

        #arr = 
        #display(Nfft2)
        arr = -1 .* collect(1:ntor) .+ (Nfft2 + 1)
        arr2 =  collect(ntor:-1:1) .+ 1
        #println(arr)
        #println(arr2)
        idxlist = [[1] ; arr ; arr2]
        #display(idxlist)
        
        #display(fftcos[1:mpol+1, idxlist])

        cn[1:mpol+1, :] = fftcos[1:mpol+1, idxlist]
        sn[1:mpol+1, :] = fftsin[1:mpol+1, idxlist]
    elseif mpol >= mpol_data && ntor < ntor_data
        arr = -1 .* collect(1:ntor) .+ (Nfft2 + 1)
        arr2 =  collect(ntor:-1:1) .+ 1
        idxlist = [[1] ; arr ; arr2]

        cn[1:mpol_data+1, :] = fftcos[1:mpol_data+1, idxlist]
        sn[1:mpol_data+1, :] = fftsin[1:mpol_data+1, idxlist]
    else
        display("rfft2d no work.")
    end
    #display(cn)
    return cn, sn


end


function irfft1D(cos_in::Array{Float64, 2}, sin_in::Array{Float64, 2}, nfft_multiplier=2)

    dim1, dim2 = size(cos_in)
    Nfft = nfft_multiplier * (dim2 - 1)
    #display(dim1)
    #display(dim2)
    #display(Nfft * 2)
    ffft = (cos_in .- 1im .* sin_in) .* Nfft

    
    #may need to change this
    #this is differretn, but I think because this is 1d it is ok.
    #will be a problemo later.
    #also do not understand why this is occuring.
    ffft[:, 1] .*= 2

    #println(ffft)
    #probbably need to pick the dim properly. Not sure what shape cos_in and sin_in have

    #need to pad ffft with zeros as it is not long enough for the inverse fft of the desired length
    #python does this automatically!
    Npad = floor(Int64,  Nfft / 2 - dim2/2) + 1
    #ft_pad = [zeros(Npad); ffft; zeros(Npad)];

    #display(Npad)

    #ok so this matches python, bit odd tbh.
    #so 2*Npad ; ffft with ifftshift makes this match python
    #I am sure there are other combos that would work.
    #so does ffft; 2*Npad, without anyshift.
    #pretty hekin annoying. But atleast we are there now!
    #display(size(ffft))
    ft_pad = zeros(ComplexF64, (dim1, dim2 + 2 * Npad))# + 2 * Npad))
    #ft_pad = [ffft; zeros(dim1, Npad); zeros(dim1, Npad)];

    ft_pad[1:end, 1:dim2] = ffft
    #may need to shift this!
    #display(length(ft_pad))
    #println(irfft(fftshift(ft_pad), 2 * Nfft) .* (Nfft / length(ffft)))
    #println(irfft(ifftshift(ft_pad), 2 * Nfft))# .* (Nfft / length(ffft)))
    #println(irfft(ft_pad, 2 * Nfft))# .* (Nfft / length
    #perhaps irfft is different in other ways???
    return irfft(ft_pad, 2 * Nfft, [2])
    #return irfft(ffft, 2 * Nfft, [2])

end

