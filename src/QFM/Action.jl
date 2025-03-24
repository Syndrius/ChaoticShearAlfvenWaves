
#copied from earlier, complete disaster
#not integrated yet.
using NLsolve #unsure if this is the best method
using FFTW

function action(pp, qq, prob::ProblemT)

    #pp is poloidal orbit periodicty, what we would describe as m
    #qq is toroidal orbit periodicty, aka n

    #pp and qq are terible names form this.

    #may need to work with iota rather than q.

    #doesn't really make sense but oh well.
    sguess = 0.6 #may need to be an arg.
    
    

    #MM = 4 # no idea
    #pqNtor = 8 #no idea
    MM = 4 #this has to be 2 * nfft_multiplier, which is defaulted to 2 in rfft functions.
    #but this at 4 and pqNtor the function is v quick.
    pqNtor = 8
    pqMpol = 24 #no idea
    #pqNtor = 2
    #pqMpol = 4 #no idea
    iota = pp/qq #iota not q!

    qN = qq * pqNtor # looks to be the number toroidal points to consider??? i.e. with a f_quad factor.
    fM = MM * pqNtor  #unsure what this is, seems to be an fft_multiplier, but scaled multiple times??

    qfM = qq * MM * pqNtor #no idea

    Nfft = MM * qq * pqNtor #no idea, guess number of points per mode by number of modes?

    #distance between toroidal points?
    dz = 2 * π / (MM * pqNtor)

    #theta distance bwteen action curves???
    dt = (2 * π  / qq) / fM

    #display(dz)
    #display(dt)

    #qN + 1 in og python
    nlist = range(0, qN)
    #don't think zeta is the toroidal component.
    ζ = range(0, Nfft-1) .* dz

    nzq = nlist .* ζ' ./ qq

    cnzq = cos.(nzq)
    snzq = sin.(nzq)

    #not certain we will actually need both sin and cos for both, seems like that is only when there is less symmetry???

    #we will almost certainly want to change the shape of these things.
    rcosarr = zeros(fM, qN+1)
    #t for theta I think, although there are multiple theta's going around.
    tcosarr = zeros(fM, qN+1)

    rsinarr = zeros(fM, qN+1)

    tsinarr = zeros(fM, qN+1)

    #absolutely no idea what this is doing...
    nvarr = zeros(fM)

    #for now we will use the default jacobian behaviour...
    #translating pythons root finding to julia may be a bit tricky regarding the jacobian part.
    #think this is just a matter of using autodiff or not...
    #if needed we can use SciPy.jl which can call scipy from julia, then we can use the same function.
    #we may need to use IntervalRootFinding.jl because we are not using a scalar function.


    #dθ = (2*π / qq) /fM #this is the theta spacing between each of the action curves to be considered
    #exactly what action curves are... not sure.

    for jpq in 0:fM-1 #not really sure what we are looping over here.

        a = jpq * dt #unsure.

        if jpq == 0
            nv0 = 0
            rcos0 = zeros(qN+1)
            tcos0 = zeros(qN+1)
            rsin0 = zeros(qN+1)
            tsin0 = zeros(qN+1)

            rcos0[1] = sguess
            tcos0[1] = 0
        else 
            #these may need to be copied!
            nv0 = nvarr[jpq]
            rcos0 = rcosarr[jpq, :]
            tcos0 = tcosarr[jpq, :]
            rsin0 = rsinarr[jpq, :]
            tsin0 = tsinarr[jpq, :]

            #rcos0[1] = sguess
            tcos0[1]  += dt
        end

        xx0 = pack_dof(nv0, rcos0, tsin0, rsin0, tcos0)

        #xx0[2] = 2.0 + 1 #+1 for packing!
        #display(xx0[1:5])
        #xx0[end] = 8.0

        #xx0[3] = 4.0
        #xx0[4] = 5.0

        #ftest = action_gradient(xx0, pp, qq, a, qN, ζ, nlist)

        #println(ftest)

        #display(ff)
       

        #short hand so it is a function of one thing for the solve
        a!(ff, xx) = action_gradient(ff, xx, pp, qq, a, qN, ζ, nlist, prob)
        #a!(ff, xx) = action_grad!(ff, xx, pp, qq, a, qN, ζ, nlist, prob)

        sol = nlsolve(a!, xx0) #etc
        #ff = ag(xx0)
        #display(ff)#[2:6])

        

        #probably want to check for convergence of solution here first!
        nv, rcos, tsin, rsin, tcos = unpack_dof(sol.zero, qN)

        #display(rcos)

        nvarr[jpq+1] = nv #this seems to be completly pointless.
        rcosarr[jpq+1, :] = rcos
        tsinarr[jpq+1, :] = tsin
        rsinarr[jpq+1, :] = rsin
        tcosarr[jpq+1, :] = tcos
    end

    r = irfft1D(rcosarr, rsinarr)
    z = LinRange(0, 2*qq*π, size(r)[end]+1)[1:end-1]



    tcosarr[:, 1] .= 0.0

    

    t = irfft1D(tcosarr, tsinarr)

    #display(t)

    r2D_alpha = zeros((qfM, Nfft))
    t2D_alpha = zeros((qfM, Nfft))

    #this indexing is cooked af.
    for i in 0:qq-1
        idx = mod(pp*i, qq)

        #display(i)

        #display(r[:, 1: i * pqNtor * MM])
        #display(r[:, 1+ i * pqNtor * MM : end])
        #display(r[:, 1: i * pqNtor * MM])
        

        r2D_alpha[1+idx * fM : (idx+1)*fM, 1 :(qq-i) * pqNtor * MM] = r[:, 1+ i * pqNtor * MM : end]

        r2D_alpha[1+idx * fM : (idx+1)*fM, 1 + (qq-i) * pqNtor * MM : end] = r[:, 1: i * pqNtor * MM]

        t2D_alpha[1+idx * fM : (idx+1)*fM, 1 :(qq-i) * pqNtor * MM] = t[:, 1+i * pqNtor * MM : end]

        t2D_alpha[1+idx * fM : (idx+1)*fM, 1 + (qq-i) * pqNtor * MM : end] = t[:, 1: i * pqNtor * MM]
    
    end

    #display(r2D_alpha) 
    #display(t2D_alpha)

    r2D_vartheta = zeros((qfM, MM * pqNtor))
    t2D_vartheta = zeros((qfM, MM * pqNtor))

    for i in 0:MM * pqNtor-1
        #v odd that this is required. but otherwise the mod function removes some of the input???
        arr = 0:qfM
        idx = @. mod(arr - i * pp, qfM) .+ 1

        


        idx = idx[1:end-1]
        #println(idx)

        #this will not work.
        r2D_vartheta[:, i+1] = r2D_alpha[idx, i+1]
        t2D_vartheta[:, i+1] = t2D_alpha[idx, i+1]
    end

    #display(t2D_vartheta)
    #odd to change r to s here but leave t as theta
    scos_surf, ssin_surf = rfft2D(r2D_vartheta, pqMpol, pqNtor)
    tcos_surf, tsin_surf = rfft2D(t2D_vartheta, pqMpol, pqNtor)

    return scos_surf, tsin_surf, ssin_surf, tcos_surf

end

#perhaps we blindly implement, despite this not really making any sense....
#alternatively we could try and get PythonCall to work in Julia, however, this sounds awful...
#function actually used for root finding!
function action_gradient(δS, xx, pp, qq, a, qN, ζ, nlist, prob)

    iota = pp / qq

    nv, rcos, tsin, rsin, tcos = unpack_dof(xx, qN)

    #display(xx[2])
    #display(rcos)
    #display(length(tsin))
    #display(length(tcos))
    #last var is unspecified atm   
    r = irfft1D(rcos, rsin)#, nfft_multiplier)

    

    #this is clearly not working as it should, as only for specific inputs
    #do the arrays match.
    t = irfft1D(tcos, tsin)

    #display(tcos)
    

    #do not have this here.
    z = ζ
    #display(length(ζ))
    #display(iota)
    #display(length(t))
    t .+= iota .* z

    #display(t)
    #println(t)

    #v unclear what this is.
    #perhaps area is fixed over the change and this is the og area?
    area = tcos[1]

    #this function is called over and over, this is terrible
    
    Br = zeros(length(r))
    Bt = zeros(length(t))
    Bz = zeros(length(z))
    met = MetT()
    B = BFieldT()

    #r, t, z are all 1d arrays. hopefully this continuos to be true.
    #this is awful.
    #unsure if there is any other way around this though!
    #but given this is computed heaps of times this will suck.
    #so looks like this is actually the problemo,
    #not the flux vs r debate.
    for i in 1:1:length(r)
        #perhaps r, t, z are not what we think they are? 
        #pretty hekin likely I think
        #nah cause the test compute B uses them as is?
        prob.compute_met(met, r[i], t[i], z[i], prob.geo.R0)

        compute_B!(B, met, prob.q, prob.isl, prob.isl2, r[i], t[i], z[i])

        #unsure how jacobian sits with all of this.
        Br[i] = B.B[1]
        Bt[i] = B.B[2]
        Bz[i] = B.B[3]


    end
    #display(Br)
    #exit()
    #display(met.gl)


    #will need to implement this!!!! -> eventually this will be linked to our og problemo.
    #may also need to scale this by the jacobian, I think we will based on our other work.
    #although we will probably start in a slab, with J = 1
    #fk this will be a nightmare.
    #Brn, Bt, Bz = test_compute_B(r, t, z)


    #println(Br[1:5])
    #println(Brn[1:5])
    #println((Br .- Brn)[1:5])
    #return

    rhs_tdot = @. Bt / Bz
    
    #nv may be the biggest mystery atm.
    #it may be the lagrange multiplier that is mentioned?
    rhs_rdot = @. Br / Bz - nv / Bz

   
    #println(rhs_rdot)
    #println(rhs_tdot)

    rhs_rdot_fft_cos, rhs_rdot_fft_sin = rfft1D(rhs_rdot)
    rhs_tdot_fft_cos, rhs_tdot_fft_sin = rfft1D(rhs_tdot)

    #println(rhs_tdot_fft_cos)
    #display(rhs_rdot)

    #xx should still be 1d. may need to be careful, think it will change once we implement the gradient.
    #ff = zeros(length(xx))

    #display(length(ff))
    #display(length(tcos))
    #pack dof does this.
    #[[nv]; rcos ; tsin[2:end] ; rsin[2:end] ; tcos] .+ 1

    δS[1] = area - a
    #don't have nlist yet!
    #display(length(rcos))
    #display(length(nlist))
    #display(length(rhs_rdot_fft_sin[1: qN + 1]))
    #println(rhs_rdot_fft_cos[1:qN+1])
    #println(nlist)
    δS[2 : qN + 2] = @. rsin * nlist / qq - rhs_rdot_fft_cos[1: qN + 1]

    δS[qN+3: 2*qN + 2] = @. (-rcos * nlist / qq - rhs_rdot_fft_sin[1:qN+1])[2:end]

    δS[2*qN+3 : 3*qN + 3] = @. (tsin * nlist / qq - rhs_tdot_fft_cos[1 : qN + 1])

    δS[2 * qN + 3] += iota #wot.

    δS[3 * qN + 4 : end] = @. (-tcos * nlist / qq - rhs_tdot_fft_sin[1:qN+1])[2:end]

    #display(ff[2:6])
    #ff looks to be ok.
    #println(ff[3 * qN + 4 : end])
    #println( @. rsin * nlist / qq - rhs_rdot_fft_cos[1: qN + 1])
    #println(ff)

    #return ff
end




function pack_dof(nv, rcos, tsin, rsin, tcos)

    #nv is just a single value I think.

    #no idea why + 1 is needed.
    #we have changed this for the 2d case.
    return [[nv]; rcos ; tsin[2:end] ; rsin[2:end] ; tcos] .+ 1

end


function pack_dof2D(nv, rcos, tsin, rsin, tcos)

    #nv is just a single value I think.

    #no idea why + 1 is needed.
    #we have changed this for the 2d case.
    return [nv; rcos ; tsin[2:end, :] ; rsin[2:end, :] ; tcos] .+ 1
    #return [nv rcos  tsin[2:end, :]  rsin[2:end, :]  tcos] .+ 1

end

function unpack_dof(xx, qN)



    nv = xx[1] - 1
    rcos = xx[2:qN + 2] .- 1

    tsin = [[0] ; xx[qN + 3:2 * qN + 2] .- 1]
    rsin = [[0] ; xx[2*qN + 3 : 3 * qN + 2] .- 1]
    tcos = xx[3 * qN + 3 : end] .- 1
    
    return nv, rcos, tsin, rsin, tcos
end

#copying, may be possible to do this more directly

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

function rfft1D(f)


    Nfft = size(f)[end]
    ffft = rfft(f)

    #println(ffft)

    cos_out = real.(ffft) ./ Nfft .* 2

    
    #1d case
    cos_out[1] /= 2
    

    sin_out = @. -imag(ffft) / Nfft * 2

    sin_out[1] = 0

    #print(sin_out)
    return cos_out, sin_out


end


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

function irfft2D(cos_in, sin_in, nfft_θ, nfft_ζ)

    mpol = size(cos_in)[1] - 1
    ntor = (size(sin_in)[2]-1) ÷ 2

    #tells us that we will need to test this with a non-zero cos_in.
    #display(cos_in)
    #display(mpol)
    #display(ntor)

    mpol_new = nfft_θ ÷ 2
    ntor_new = nfft_ζ ÷ 2

    cos_pad = zeros(mpol_new+1, 2*ntor_new)
    sin_pad = zeros(mpol_new+1, 2*ntor_new)

    #display(-1 .* collect(1:ntor) .+1)
    #display(collect(ntor:0))
    arr = -1 .* collect(1:ntor) .+ size(cos_in)[2] .+ 1
    arr2 =  collect(ntor:-1:1) .+ 1
    #println(arr)
    #println(arr2)
    #this is wrong for sure.
    idxlist = [[1] ; arr ; arr2]
    #display(idxlist)

    cos_pad[1:mpol+1, idxlist[:]] = cos_in
    sin_pad[1:mpol+1, idxlist] = sin_in

    cos_pad[1, :] .*= 2
    cos_pad[end, :] .*= 2
    sin_pad[1, :] *= 2
    sin_pad[end, :] .*= 2

    #display(cos_pad)

    #so this works for some reason...
    #need to figure this shit out.
    #so d needs to be the new size of the array.
    #v unclear what is deciding that!
    #not sure if this works in general! #maybe it does
    #also not sure if we need to do the fft with dimensions flipped.
    fout = irfft(cos_pad .- 1im .* sin_pad, 2*mpol_new, [1, 2])

    return fout .* nfft_θ .* nfft_ζ ./ 2

end


#function used as a test for comparison with Zhisong
#want to replace with our normal B method.
function test_compute_B(r, t, z)

    #q = @. 2 / r^2

    #zhisongs names for this are frankly unacceptable, we have 
    #changed to normal coord names.
    k = 0.0018 #same as python example notebook
    #gotta deal with vectors here big rip.
    #we will need to make this conform to our normal method later!
    #unsure if this will always have vector inputs.
    #Br = zeros(length(r), length(t), length(z))
    #Bt = zeros(length(r), length(t), length(z))
    #Bz = zeros(length(r), length(t), length(z))
    #this will cause some problemos later, 
    #this coord setup is weird,
    #Br is a 1d array.
    Br = @. - k * (sin(2 * t - z) + sin(3 * t - 2 * z))
    Bt = @. r #/ q
    Bz = ones(length(r)) #r

    return Br, Bt, Bz
end
