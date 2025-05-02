"""
So this file is starting to come together.
Big issues still:
- What is field wrapping? Why are we doing it?
- Parameters/shape of arrays, need to understand the Rfft stuff more.
- The actual equations we are solving.
- Doubly so for JM, that makes little to no sense.
- JM is still allocating a lot. Understanding the above will help with that.
"""
#we have basic functionality again. This is all still pretty cooked though!

#no idea what to call this damn thing!
#not even sure we want this tbh!
#note we are incosisntently using coefs and coeffs
struct CoefficientsT
    nv :: Array{Float64, 1} #this is actually just a float. but we want to mutate it.
    rcos :: Array{Float64, 1}
    θcos :: Array{Float64, 1}
    rsin :: Array{Float64, 1}
    θsin :: Array{Float64, 1}
    function CoefficientsT(N::Int64)
        new(zeros(1), zeros(N), zeros(N), zeros(N), zeros(N))
    end
end


#struct to condense the inputs for action grad. allowing preallocation of variables.
#perhaps this should be split up into multiple structs. this is wild.
#the split between this struct and the next is much more arbitrary than we would like.
struct ActionGradInputsT
    p :: Int64
    q :: Int64
    Ntor :: Int64
    r :: Array{Float64}
    θ :: Array{Float64}
    ζ :: StepRangeLen{Float64}
    prob :: ProblemT
    met :: MetT
    B :: BFieldT
    nlist :: Array{Int64}
    rdot :: Array{Float64}
    θdot :: Array{Float64}
    function ActionGradInputsT(p::Int64, q::Int64, Ntor::Int64, fft_size::Int64, ζ::StepRangeLen{Float64}, prob::ProblemT, met::MetT, B::BFieldT, nlist::UnitRange{Int64})
        #perhaps the inputs are a bit restrictive...
        new(p, q, Ntor, zeros(2*fft_size), zeros(2*fft_size), ζ, prob, met, B, nlist, zeros(2*fft_size), zeros(2*fft_size))
    end
end

#struct for storing the data needed for efficient fourier transform
struct FTDataT
    ift_1D_p :: AbstractFFTs.ScaledPlan
    rft_1D_p :: FFTW.rFFTWPlan
    ift_r1D :: Array{ComplexF64}
    ift_θ1D :: Array{ComplexF64}
    rdot_fft_cos :: Array{Float64}
    rdot_fft_sin :: Array{Float64}
    θdot_fft_cos :: Array{Float64}
    θdot_fft_sin :: Array{Float64}
    ft_r1D :: Array{ComplexF64}
    ft_θ1D :: Array{ComplexF64}
    function FTDataT(fft_size::Int64)
        #obvs a bit dumb that we are adding 1 to everything!
        temp = zeros(2*fft_size)
        rft_1D_p = plan_rfft(temp)
        temp = zeros(ComplexF64, fft_size+1)
        irft_1D_p = plan_irfft(temp, 2*fft_size)
        new(irft_1D_p, rft_1D_p, zeros(ComplexF64, fft_size+1), zeros(ComplexF64, fft_size+1), zeros(fft_size+1), zeros(fft_size+1), zeros(fft_size+1), zeros(fft_size+1), zeros(ComplexF64, fft_size+1), zeros(ComplexF64, fft_size+1))
    end
end


#struct to condense the inputs for action grad jm. allowing preallocation of variables.
#this takes exactly the same args. Don't think it should!
#still needs work, we will need to fix action_grad_jm first tbh.
#this is even more of a disaster. Probably want the normal one and an extra one, just so this isn't quite so fkn cooked.
struct ActionGradInputsJMT
    ooBζ :: Array{Float64}
    ooBζcos :: Array{Float64}
    ooBζsin :: Array{Float64}
    cosnzq :: Array{Float64, 2}
    sinnzq :: Array{Float64, 2}
    rdot_dr :: Array{Float64}
    rdot_dθ :: Array{Float64}
    θdot_dr :: Array{Float64}
    θdot_dθ :: Array{Float64}
    rdot_drcos :: Array{Float64, 2}
    rdot_drsin :: Array{Float64, 2}
    rdot_dθcos :: Array{Float64, 2}
    rdot_dθsin :: Array{Float64, 2}
    θdot_drcos :: Array{Float64, 2}
    θdot_drsin :: Array{Float64, 2}
    θdot_dθcos :: Array{Float64, 2}
    θdot_dθsin :: Array{Float64, 2}
    #rdot :: Array{Float64, 2} #no fkn idea what the shape of this is tbh!
    #θdot :: Array{Float64, 2} #no fkn idea what the shape of this is tbh!
    function ActionGradInputsJMT(fft_size::Int64, ip::ActionGradInputsT)
        nzq = ip.nlist .* ip.ζ' ./ ip.q
        M = length(ip.nlist)
        N = 2 * fft_size

        #fft_size is probably not the best name for this.
        #note that we don't actually know the size of ooBζcos!
        new(zeros(N), zeros(fft_size+1), zeros(fft_size+1), cos.(nzq), sin.(nzq), zeros(N), zeros(N), zeros(N), zeros(N), zeros(M, N), zeros(M, N), zeros(M, N), zeros(M, N), zeros(M, N), zeros(M, N), zeros(M, N), zeros(M, N))
    end
end


#not sure abot this name tbh. It is more than just the action.
#may be worth splitting, as the start is more naturally called the action, then we can have a 'wrapping' or something function.
function action(p::Int64, q::Int64, prob::ProblemT, met::MetT, B::BFieldT, MM::Int64, M::Int64, N::Int64, sguess=0.5::Float64)#, sguess::Float64, MM::Int64, M::Int64, N::Int64)

    #we take our r, θ original toroidal coordinates as an arbitrary fourier expansion in terms of zeta to denote our surface
    #r(ζ) = r_0^c ∑_{n=1}^qN (r_n^c cos(nζ/q) + r_n^s sin(nζ/q))
    #θ(ζ) = θ_0^c ∑_{n=1}^qN (θ_n^c cos(nζ/q) + θ_n^s sin(nζ/q))

    #M is fourier resolution in poloidal direction
    #N is the fourier resolution in in toroidal direction.
    #still pretty unsure exactly what this is.
    #MM = 4 #this has to be 2 * nfft_multiplier, which is defaulted to 2 in rfft functions.

    #number of toroidal modes.
    qN = q * N
    fM = MM * N #the number of psuedo fieldlines to be found.
    #number of poloidal points once the fiedlines are wrapped. Eg this is the extension of the domain from 2π to 2qπ.
    qfM = q * fM
    #Number of toroidal points for each action gradient calculation, unsure if this is actually meant to be the same as qfM.
    Nfft= qfM

    Ntor = qN + 1 #actual size of our fft arrays.

    #no idea what these are either!
    dζ = 2π / fM
    dθ = 2π / qfM

    nlist = range(0, qN)
    #don't think zeta is the toroidal component.
    ζ = range(0, Nfft-1) .* dζ

    #array that stores the all the coefficients for root finding
    #-2 because we don't include the sin 0 coeffs.
    x0 = zeros(4*Ntor+1-2)

    #these store the minimised values for the coefficeints at each iteration
    rcosarr = zeros(fM , Ntor)
    rsinarr = zeros(fM , Ntor)
    θcosarr = zeros(fM , Ntor)
    θsinarr = zeros(fM , Ntor)
    nvarr = zeros(fM)
    
    nfft_mulitplier = 2
    #qunatity used to define the size of most arrays needed.
    size_of_fft_array = nfft_mulitplier * (Ntor - 1)
    #creates a struct with prealocated memory as input
    ip = ActionGradInputsT(p, q, Ntor, size_of_fft_array, ζ, prob, met, B, nlist)

    #creates a struct with preallocated memeory for fourier data.
    ftd = FTDataT(size_of_fft_array)

    ipjm = ActionGradInputsJMT(size_of_fft_array, ip)

    #struct storing the actual coefficients to be found.
    coefs = CoefficientsT(qN+1)

    #iterates over different poloidal areas, this creates multiple psuedo fieldlines that minamise the action.
    #these are then combined into a surface below.
    for jpq in 0:fM-1
        #maybe this is supposed to be alpha?
        #or area? fk knows.
        a = jpq * dθ
        if jpq == 0
            #initial guesses for the alg.
            coefs.rcos[1] = sguess
            #everything else is zero.
            rcos0 = zeros(qN+1)
            rsin0 = zeros(qN+1)
            θcos0 = zeros(qN+1)
            θsin0 = zeros(qN+1)

            rcos0[1] = sguess
            θcos0[1] = 0

            nv0 = 0
        else
            #otherwise use the last solution as a starting point.
            coefs.rcos .= rcosarr[jpq, :]
            coefs.rsin .= rsinarr[jpq, :]
            coefs.θcos .= θcosarr[jpq, :]
            coefs.θsin .= θsinarr[jpq, :]
            coefs.nv[1] = nvarr[jpq]
            rcos0 = rcosarr[jpq, :]
            rsin0 = rsinarr[jpq, :]
            θcos0 = θcosarr[jpq, :]
            θsin0 = θsinarr[jpq, :]
            nv0 = nvarr[jpq]

            θcos0[1]  += dθ 
            coefs.θcos[1] += dθ
        end
        #combine the values into a single array for the root findning alg's argument
        pack_coeffs!(x0, coefs, Ntor)

        #creates a standin function that fits the form needed for NLSolve.
        ag!(δS, x) = action_grad!(δS, x, a, coefs, ip, ftd)

        #creates a standin function that fits the form needed for NLSolve.
        #this function defines the gradient, improving efficiency.
        ag_JM!(JM, x) = action_grad_jm!(JM, x, a, coefs, ip, ftd, ipjm)

        #solve for the psuedo field line.
        sol = nlsolve(ag!, ag_JM!, x0)

        #unpacks the single array solution into CoefficientsT struct.
        unpack_coeffs!(sol.zero, coefs, Ntor)

        #stores this field line.
        rcosarr[jpq+1, :] = coefs.rcos
        rsinarr[jpq+1, :] = coefs.rsin
        θcosarr[jpq+1, :] = coefs.θcos
        θsinarr[jpq+1, :] = coefs.θsin
        nvarr[jpq+1] = coefs.nv[1]
        
    end

    return wrap_field_lines(rcosarr, rsinarr, θcosarr, θsinarr, MM, M, N, p, q, qfM, Nfft)
end



"""

Function that turns the Pseudo field lines into a surface. Probably the least clear function in the entire package.
"""
function wrap_field_lines(rcosarr, rsinarr, θcosarr, θsinarr, MM, M, N, p, q, qfM, Nfft)

    #think the params we have here and above would be simpler if we knew what any of them were.
    fM = MM * N 

    #converts the fourier coefficients back to the normal values.
    r = irfft1D(rcosarr, rsinarr)
    ζ = LinRange(0, 2*q*π, size(r)[end]+1)[1:end-1]

    θcosarr[:, 1] .= 0.0
    θ = irfft1D(θcosarr, θsinarr)


    #this defines r(α, ζ), extending the domain from 2π to 2qπ.
    r2D_alpha = zeros((qfM, Nfft))
    θ2D_alpha = zeros((qfM, Nfft))

    #this just extends the solution from 2π to 2qπ, as this is the periodicity of the pesudo fieldline.
    for i in 0:q-1
        idx = mod(p*i, q)

        r2D_alpha[1+idx * fM : (idx+1)*fM, 1 :(q-i) * N * MM] = r[:, 1+ i * N * MM : end]

        r2D_alpha[1+idx * fM : (idx+1)*fM, 1 + (q-i) * N * MM : end] = r[:, 1: i * N * MM]

        θ2D_alpha[1+idx * fM : (idx+1)*fM, 1 :(q-i) * N * MM] = θ[:, 1+i * N * MM : end]

        θ2D_alpha[1+idx * fM : (idx+1)*fM, 1 + (q-i) * N * MM : end] = θ[:, 1: i * N * MM]
    
    end

    r2D_vartheta = zeros((qfM, MM * N))
    θ2D_vartheta = zeros((qfM, MM * N))

    #here we actually get r(ϑ, φ), creating the result we want.
    for i in 0:MM * N-1
        #v odd that this is required. but otherwise the mod function removes some of the input???
        arr = 0:qfM
        idx = @. mod(arr - i * p, qfM) .+ 1

        idx = idx[1:end-1]

        r2D_vartheta[:, i+1] = r2D_alpha[idx, i+1]
        θ2D_vartheta[:, i+1] = θ2D_alpha[idx, i+1]
    end

    #converts the final results back to fourier coefficients, I guess for easier interpolation.
    scos_surf, ssin_surf = rfft2D(r2D_vartheta, M, N)
    tcos_surf, tsin_surf = rfft2D(θ2D_vartheta, M, N)

    return scos_surf, tsin_surf, ssin_surf, tcos_surf

end



"""
    action_grad!(δS::Array{Float64}, x::Array{Float64}, p::Int64, q::Int64, a::Float64, qN::Int64, ζ, nlist, prob::ProblemT)

Computes the gradient of the action for pseudo fieldlines.
"""
function action_grad!(δS::Array{Float64}, x::Array{Float64}, a::Float64, coefs::CoefficientsT, ip::ActionGradInputsT, ftd::FTDataT)

    unpack_coeffs!(x, coefs, ip.Ntor)

    iota = ip.p/ip.q

    area = coefs.θcos[1]

    get_r_t!(ip.r, ip.θ, coefs, ftd.ift_1D_p, ftd.ift_r1D, ftd.ift_θ1D, ip.Ntor)

    ip.θ .+= iota .* ip.ζ

    for i in 1:1:length(ip.r)
        #gross amount of nested structs.
        ip.prob.met(ip.met, ip.r[i], ip.θ[i], ip.ζ[i], ip.prob.geo.R0)

        #alt smaller function, not being used atm
        #compute_B!(B.B, met.J, prob.q, prob.isl, prob.isl2, r[i], θ[i], ζ[i])
        compute_B!(ip.B, ip.met, ip.prob.q, ip.prob.isls, ip.r[i], ip.θ[i], ip.ζ[i])

        #unsure if perhaps the jacobian is required in here somewhere.
        ip.θdot[i] = ip.B.B[2] / ip.B.B[3]
        #pretty sure the Jacobian is requried here, this will mess up the JM function though!
        ip.rdot[i] = (ip.B.B[1] - coefs.nv[1]) / ip.B.B[3]
    end
    

    rfft1D!(ftd.rdot_fft_cos, ftd.rdot_fft_sin, ip.rdot, ftd.ft_r1D, ftd.rft_1D_p)
    rfft1D!(ftd.θdot_fft_cos, ftd.θdot_fft_sin, ip.θdot, ftd.ft_θ1D, ftd.rft_1D_p)

    
    #note that thes equations are not clear from the paper.
    #I think these basically have nothing to do with the equations as written in the paper.

    #we are essentially saying that our current guess of the path is given by 
    #r(ζ) = r_0^c ∑_{n=1}^qN (r_n^c cos(nζ/q) + r_n^s sin(nζ/q))
    #this should have a plus ζ term though!
    #θ(ζ) = θ_0^c ∑_{n=1}^qN (θ_n^c cos(nζ/q) + θ_n^s sin(nζ/q))

    #we then take the difference between the ideal dynamics
    #∂r/∂ζ = B^r/B^ζ - ν
    #∂θ/∂ζ = B^θ/B^ζ

    #where ν is the extra area bit that allows approx solutions to be found, true, integrable dynamics would satisfy these two equation without needing the ν
    #we then find the coeffs that define our initial guess
    #based on the difference between the current guess and the `ideal' dynamics based on the field.
    #the equations quoted seem unhelpful.

    #og packing
    #[[nv]; rcos ; tsin[2:end] ; rsin[2:end] ; tcos]
    qN = ip.Ntor - 1
    δS[1] = area - a

    δS[2 : qN + 2] = @. coefs.rsin * ip.nlist / ip.q - ftd.rdot_fft_cos[1: qN + 1]

    δS[qN+3: 2*qN + 2] = @. (-coefs.rcos * ip.nlist / ip.q - ftd.rdot_fft_sin[1:qN+1])[2:end]

    δS[2*qN+3 : 3*qN + 3] = @. (coefs.θsin * ip.nlist / ip.q - ftd.θdot_fft_cos[1 : qN + 1])

    δS[2 * qN + 3] += iota #wot.

    δS[3 * qN + 4 : end] = @. (-coefs.θcos * ip.nlist / ip.q - ftd.θdot_fft_sin[1:qN+1])[2:end]
    
    #new packing
    #[rcos ; rsin[2:end] ; tcos ; tsin[2:end] ; [nv]]
    #=
    δS[1:Ntor] = @. CT.rsin * nlist / q - rdot_fft_cos[1: Ntor]
    δS[Ntor+1:2*Ntor] = @. (CT.θsin * nlist / q - θdot_fft_cos[1 : Ntor])
    δS[Ntor+1] += iota
    δS[2*Ntor+1:3*Ntor-1] = @. (-CT.θcos * nlist / q - θdot_fft_sin[1:Ntor])[2:end]
    δS[3*Ntor+1-1:4*Ntor-2] = @. (-CT.rcos * nlist / q - rdot_fft_sin[1:Ntor])[2:end]
    δS[end] = area - a
    =#
end

#computes the jacobian matrix for the action grad for efficient root finding.
#fair bit of overlap with this function and the og
#may be better to combine them
#but we are probably getting to diminishing returns.
#this is still an allocating unclear mess. But the inputs have been reduced lol.
"""

#TODO
Computes the jacobian matrix of the action gradient, allows for more efficient solving.
Some allocations have been removed, however, there are still a lot, mostly with arrays that are not clear.
"""
function action_grad_jm!(JM::Array{Float64, 2}, x::Array{Float64}, a::Float64, coefs::CoefficientsT, ip::ActionGradInputsT, ftd::FTDataT, ipjm::ActionGradInputsJMT)


    iota = ip.p/ip.q


    unpack_coeffs!(x, coefs, ip.Ntor)

    area = coefs.θcos[1]

    get_r_t!(ip.r, ip.θ, coefs, ftd.ift_1D_p, ftd.ift_r1D, ftd.ift_θ1D, ip.Ntor)
    #get_r_t!(r, θ, rcos, tsin, rsin, tcos, ift_1D_p, ft_r1D, ft_θ1D, Ntor)
    ip.θ .+= iota .* ip.ζ

    
    for i in 1:1:length(ip.r)
        
        ip.prob.met(ip.met, ip.r[i], ip.θ[i], ip.ζ[i], ip.prob.geo.R0)

        #need B.dB now, so we just compute the full thing!
        compute_B!(ip.B, ip.met, ip.prob.q, ip.prob.isls, ip.r[i], ip.θ[i], ip.ζ[i])

        #unsure how jacobian sits with all of this.
        #dBrdr/Bζ - (Br - nv) / Bζ^2 * dBζdr
        ipjm.rdot_dr[i] = ip.B.dB[1, 1] / ip.B.B[3] - (ip.B.B[1] - coefs.nv[1]) / ip.B.B[3] * ip.B.dB[3, 1]
        #dBrdθ/Bζ - (Br - nv) / Bζ^2 * dBζdθ
        ipjm.rdot_dθ[i] = ip.B.dB[1, 2] / ip.B.B[3] - (ip.B.B[1] - coefs.nv[1]) / ip.B.B[3] * ip.B.dB[3, 2]
        #dBθdr/Bζ - Bθ / Bζ^2 * dBζdr
        ipjm.θdot_dr[i] = ip.B.dB[2, 1] / ip.B.B[3] - ip.B.B[2] / ip.B.B[3] * ip.B.dB[3, 1]
        #dBθdθ/Bζ - Bθ / Bζ^2 * dBζdθ
        ipjm.θdot_dθ[i] = ip.B.dB[2, 2] / ip.B.B[3] - ip.B.B[2] / ip.B.B[3] * ip.B.dB[3, 2]
        
        ipjm.ooBζ[i] = 1 / ip.B.B[3]
    end



    #need to change this!
    ipjm.ooBζcos[:], ipjm.ooBζsin[:] = rfft1D_simple(ipjm.ooBζ)

    #no idea what these even are tbh.
    #not actually sure how big these are needed to be yet, but they should be the same shape as nzq.

    M, N = size(ipjm.rdot_drcos)
    for j in 1:N, i in 1:M
        #may need to transpose these!
        #just genuinly have no idea what is going on with the fft step.
        ipjm.rdot_drcos[i, j] = ipjm.rdot_dr[j] * ipjm.cosnzq[i, j]
        ipjm.rdot_drsin[i, j] = ipjm.rdot_dr[j] * ipjm.sinnzq[i, j]
        ipjm.rdot_dθcos[i, j] = ipjm.rdot_dθ[j] * ipjm.cosnzq[i, j]
        ipjm.rdot_dθsin[i, j] = ipjm.rdot_dθ[j] * ipjm.sinnzq[i, j]
        ipjm.θdot_drcos[i, j] = ipjm.θdot_dr[j] * ipjm.cosnzq[i, j]
        ipjm.θdot_drsin[i, j] = ipjm.θdot_dr[j] * ipjm.sinnzq[i, j]
        ipjm.θdot_dθcos[i, j] = ipjm.θdot_dθ[j] * ipjm.cosnzq[i, j]
        ipjm.θdot_dθsin[i, j] = ipjm.θdot_dθ[j] * ipjm.sinnzq[i, j]
    end


    #this is the point where Zhisongs code actually makes zero sense whatso ever.
    dummy_nv = zeros(1, size(ipjm.rdot_drcos)[2])
    #like what kind of monstrosoty is this creating.
    #no idea what the shape of this is, so no preallocation yet.
    drdot = pack_dof2D(dummy_nv, ipjm.rdot_drcos, ipjm.rdot_dθsin, ipjm.rdot_drsin, ipjm.rdot_dθcos) .- 1
    dθdot = pack_dof2D(dummy_nv, ipjm.θdot_drcos, ipjm.θdot_dθsin, ipjm.θdot_drsin, ipjm.θdot_dθcos) .- 1

    #how can we possibly want whatever the output of this is.
    #like wot the fek.
    drdot_cos, drdot_sin = rfft1D_JM(drdot)
    dθdot_cos, dθdot_sin = rfft1D_JM(dθdot)

    qN = ip.Ntor - 1

    JM .= 0.0 #needed as they are not set to zero by default

    #unsure what any of this even is tbh. Not even sure what the JM is representing here.
    #now reasnably confident these are working.
    #genuuinly now way this can be right,
    #should be JM[coef, deriv]
    #these are implying that the derivative of the rcos with respect to every coef is the same?
    #probbaly an indexing issue.
    JM[2:qN+2, 1:end] .= -drdot_cos[1:qN+1, 1:end]
    JM[qN+3:2*qN+2, 1:end] .= -drdot_sin[2:qN+1, 1:end]

    JM[2*qN+3:3*qN+3, 1:end] .= -dθdot_cos[1:qN+1, 1:end]
    JM[3*qN+4:end, 1:end] .= -dθdot_sin[2:qN+1, 1:end]

    JM[1, 3*qN+3] += 1  #θcos[0] apparantly

    #unsure about this, python has swapped to using arange for some reason.
    #given up on trying to get the indexing right for now.
    #this obvs won't work either lol.
    #this needs to be tested.
    #python indexing for this bit is kind of weird
    #unsure how best to replicate it.
    #CartesianIndex here creates pairwise indicies from two lists
    #eg ([1, 2, 3], [2, 3, 4]) -> [1, 2], [2, 3], [3, 4]
    JM[CartesianIndex.(2:qN+2, 2*qN+2:3*qN+2)] += ip.nlist / ip.q #rsin

    JM[CartesianIndex.(qN+3:2*qN+2, 3:qN+2)] += -(ip.nlist / ip.q)[2:end] #-rsin[2:end]

    JM[CartesianIndex.(2*qN+3:3*qN+3, qN+2:2*qN+2)] += ip.nlist / ip.q #θsin

    JM[CartesianIndex.(3*qN+4:4*qN+3, 3*qN+4:4*qN+3)] += (-ip.nlist / ip.q)[2:end] #-θcos[2:end]

    #then derivative w.r.t nv
    #srs what even is this.
    JM[2:qN+2, 1] += ipjm.ooBζcos[1:qN+1]
    JM[qN+3:2*qN+2, 1] += ipjm.ooBζsin[2:qN+1]
    
end



"""

Function that packs the CoefficientsT struct into a 1d array for solving.
"""
function pack_coeffs!(x::Array{Float64}, coefs::CoefficientsT, N)

    #I think we will still want to change the order of this.
    #just need to be a wee bit careful.
    #should be easier with this function and the struct.
    #x[1] = CT.nv[1]
    #changed the order hopefully not a mistake!
    #x[1:N] = CT.rcos[:]
    #x[N+1:2*N-1] = CT.rsin[2:end]
    #x[2*N-1+1:3*N-1] = CT.θcos[:]
    #x[3*N-1+1:4*N-1-1] = CT.θsin[2:end]
    #x[end] = CT.nv[1]
    #so the plus 1 is extremely important for efficient solving!
    #why??>>?>
    #og packing
    x .= [coefs.nv ; coefs.rcos ; coefs.θsin[2:end] ; coefs.rsin[2:end] ; coefs.θcos] .+ 1
    #new packing
    #x .= [CT.rcos ; CT.rsin[2:end] ; CT.θcos ; CT.θsin[2:end] ; CT.nv] .+ 1

end


"""

Function that unpacks the 1d list of coefficients into the CoefficientsT struct.
"""
function unpack_coeffs!(x::Array{Float64}, coefs::CoefficientsT, N::Int64)

    #I think we will still want to change the order of this.
    #just need to be a wee bit careful.
    #should be easier with this function and the struct.
    #x[1] = CT.nv[1]
    
    #changed the order hopefully not a mistake!
    #=
    CT.rcos[:] = x[1:N] .- 1
    CT.rsin[2:end] = x[N+1:2*N-1] .- 1
    CT.θcos[:] = x[2*N-1+1:3*N-1] .-1
    CT.θsin[2:end] = x[3*N-1+1:4*N-1-1] .-1
    CT.nv[1] = x[end] - 1
    CT.rsin[1] = 0.0
    CT.θsin[1] = 0.0
    =#

    #og packing
    
    qN = N - 1
    coefs.nv[1] = x[1] - 1
    coefs.rcos .= x[2:qN + 2]  .- 1

    coefs.θsin .= [[0] ; x[qN + 3:2 * qN + 2] .- 1]
    coefs.rsin .= [[0] ; x[2*qN + 3 : 3 * qN + 2] .- 1]
    coefs.θcos .= x[3 * qN + 3 : end] .- 1
    
end

#this needs to be changed/removed
#so this is still being used!
#hard to fix this when we don't know why this is being done, or why the result is fourier transformed. Seems like complete madness.
function pack_dof2D(nv, rcos, tsin, rsin, tcos)

    #nv is just a single value I think.

    #no idea why + 1 is needed.
    #we have changed this for the 2d case.
    return [nv; rcos ; tsin[2:end, :] ; rsin[2:end, :] ; tcos] .+ 1
    #return [nv rcos  tsin[2:end, :]  rsin[2:end, :]  tcos] .+ 1

end



########################################################

#function for simplified computing of B, easy comparison with Zhisong's python code.
function test_compute_B!(Br, Bθ, Bζ, r, t, z)

    q = @. 2 / r^2

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
    @. Br = - k * (sin(2 * t - z) + sin(3 * t - 2 * z))
    @. Bθ = r #/ q
    @. Bζ = r #ones(length(r))
    Bζ .= 1.0

end
