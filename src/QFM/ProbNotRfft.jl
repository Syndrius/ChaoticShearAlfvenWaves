"""
    rfft1D!(cos_coef::Array{Float64}, sin_coef::Array{Float64}, func::Array{Float64}, fft_res::Array{ComplexF64}, plan::FFTW.rFFTWPlan, Nfft::Int64)

Takes the real fft in 1D in place. fft_res is a temporary array to store the rfft, before it is converted sin and cos.
For func of size N, fft_res, cos_coef and sin_coef must all be size N/2 +1, while the plan is created based on size N array.
"""
function rfft1D!(cos_coef::Array{Float64, 1}, sin_coef::Array{Float64, 1}, func::Array{Float64, 1}, fft_res::Array{ComplexF64, 1}, plan::FFTW.rFFTWPlan, Nfft::Int64)


    #takes the rfft.
    fft_res .= plan * func

    #note the additional scaling by a*π here has been divided in the main equations.

    #gets the cos component, 
    cos_coef .= real.(fft_res) ./ Nfft .* 2
    #this is required to reflect the *= 2 in irfft1D.
    cos_coef[1] /= 2

    #gets the sin component
    sin_coef .= @. -imag(fft_res) / Nfft * 2
    #sin(n=0) term is always zero
    sin_coef[1] = 0.0

end

"""

Computes the 1d rfft for a 2 dimensional input array. Used in action_grad_j!, when coefficients are combined into a single array for efficiency.
"""
function rfft1D!(cos_coef::Array{Float64, 2}, sin_coef::Array{Float64, 2}, func::Array{Float64, 2}, fft_res::Array{ComplexF64, 2}, plan::FFTW.rFFTWPlan, Nfft::Int64)


    #takes the rfft.
    fft_res .= plan * func

    #note the additional scaling by a*π here has been divided in the main equations.

    #gets the cos component, 
    cos_coef .= real.(fft_res) ./ Nfft .* 2
    #this is required to reflect the *= 2 in irfft1D.
    cos_coef[1, :] ./= 2

    #gets the sin component
    sin_coef .= @. -imag(fft_res) / Nfft * 2
    #sin(n=0) term is always zero
    sin_coef[1, :] .= 0.0

end

"""
    irfft1D!(res::Array{Float64}, cos_in::Array{Float64}, sin_in::Array{Float64}, irfft_p::AbstractFFTs.ScaledPlan, fft_arr::Array{ComplexF64}, qN::Int64, nfft::Int64)

Takes the inverse real fourier transform. cos_in and sin_in are length qN+1, temporary array to be irfft'd is length qN+1, and res is size 2*qN*nfft, irfft plan is based on size of res.
"""
function irfft1D!(res::Array{Float64}, cos_in::Array{Float64}, sin_in::Array{Float64}, irfft_p::AbstractFFTs.ScaledPlan, fft_arr::Array{ComplexF64}, qN::Int64, nfft::Int64)

    #scaled size of array
    Nfft = nfft * qN

    fft_arr .= 0.0

    #reconstruct the fourier representation
    #fft_arr[1:qN+1] .= @. (cos_in - 1im * sin_in) * Nfft
    fft_arr[1:qN+1] .= @. (cos_in - 1im * sin_in) * Nfft

    #don't think this is needed either!
    #so for some reason this is essential, probably explains the /2 for cos_coeff[1]
    fft_arr[1] *= 2

    res .= irfft_p * fft_arr
end


#used by wrap_field_lines, unclear exactly what is happening
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
    #display(size(ft_pad))
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


#used by wrap_field_lines, bit unclear
function rfft2D(f, mpol, ntor)

    #zhisongs code considers case where mpol and ntor are unspecified.
    #variable M is only used here, and it seems like it won't work half the time anyway
    #think it should be removed tbh!
    #perhaps it is just specifiying how many poloidal points we want in our final surfaces?
    #seems like that should be dictated by something other than a random af input?

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
