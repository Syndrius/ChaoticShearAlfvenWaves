"""

This file and GradAction.jl extremize the action
∮ A⋅dl - ν(1/2aπ ∮θ∇ζ⋅dl - bπ - α)
to find qfm surfaces.

This involves finding the zeros of the gradient of the action, written as 3 equations,
1) ṙ - B^r/B^ζ + ν/2aπJB^ζ = 0
2) θ̇ - B^θ/B^ζ = 0
3) 1/2aπ∫_0^2aπ θ dζ - bπ - α = 0

This is done by expanding r, θ via their fourier coefficients,

r(ζ) = r_0^c + ∑_{n=1}^aN (r_n^c cos(nζ/q) + r_n^s sin(nζ/q))
θ(ζ) = θ_0^c + ι*ζ + ∑_{n=1}^aN (θ_n^c cos(nζ/q) + θ_n^s sin(nζ/q))

leading to 4*aN+3 eqautions to solve (+3 due to n=0 terms and ν equation).

This is done with the NLSolve package.
Following their notation, in GradAction.jl we have the functions
action_grad_f! - Computes the residuals of the above equations.
action_grad_j! - Computes the jacobian of the above equations.
action_grad_fj! - Computes both, used for efficient solving.

Each function required many inputs which have been preallocated for efficieny.
These inputs are stored in the Structs below, which also follow the naming pattern used by NLSolve.

We use the terminology `rhs' for the parts of the equations above that are not ṙ or θ̇, eg ṙ_rhs is B^r/B^ζ - ν/2aπJB^ζ.

The equations are stored in the array δS, the first aN+1 elements are from the first equation above, multiplied by cos and integrated, to extract a single coeffienct, in this case the r_n^s coefficeints due to the derivative.

The next aN equations are the same equation but with sin.
Then aN+1 from the second equation with cos followed by aN equations from sin.
Finally, the remaining equation 3).

The Jacobian matrix used requires that derivatives of each of these equations is taken with respect to all the variables.

The order of this matrix depends on the order of the variables, which we take to be
[r^c_n, r^s_n, θ^c_n, θ^s_n, ν]

Final note on the size of the arrays used,
- the coefficients are length aN (+1 for the ^c cases from n=0)
- Arrays are expanded with nfft to reduce aliasing issues, so N = nfft * aN
- These are irfft'd to get values for r,θ. these arrays are size 2*(N-1)
- This is the size of most arrays used in the computation.

"""


"""
Structure that stores the variables to solve for.
This includes the fourier coefficients for r and θ, and the lagrange multiplier ν.
"""
struct CoefficientsT
    rcos :: Array{Float64, 1}
    θcos :: Array{Float64, 1}
    rsin :: Array{Float64, 1}
    θsin :: Array{Float64, 1}
    ν :: Array{Float64, 1} #this is actually just a float. but we want to mutate it.
    function CoefficientsT(N::Int64)
        new(zeros(N), zeros(N), zeros(N), zeros(N), zeros(1))
    end
end


"""
Structure for storing data used for the fourier transformations.
This includes temp arrays to store the output of real fft's before being split into sin/cos components and fourier transform plans.
"""
struct FTDataT
    irfft_1D_p :: AbstractFFTs.ScaledPlan
    rfft_1D_p :: FFTW.rFFTWPlan
    temp :: Array{ComplexF64}
    jm_temp :: Array{ComplexF64, 2} #probably 2 dims.
    jmrfft_1D_p :: FFTW.rFFTWPlan 
    function FTDataT(aN::Int64, nfft::Int64)
        #size of arrays in fourier space
        fft_size = 2 * aN * nfft
        temp = zeros(fft_size)
        rfft_1D_p = plan_rfft(temp)
        #size of arrays after they have been irfft'd.
        ifft_size = aN * nfft + 1
        temp = zeros(ComplexF64, ifft_size)
        irfft_1D_p = plan_irfft(temp, fft_size)

        #size of combined array for efficient fourier transform
        comb_size = 4 * aN + 2
        temp = zeros(fft_size, comb_size)
        jmrfft_1D_p = plan_rfft(temp, [1]) 
        new(irfft_1D_p, rfft_1D_p, zeros(ifft_size), zeros(ifft_size, comb_size), jmrfft_1D_p)
    end
end

"""
Struct storing the 'rhs' data in each equation, includes the computed value and cos/sin coefficients once fft'd.
"""
struct RHST
    rhs :: Array{Float64}
    cn :: Array{Float64}
    sn :: Array{Float64}
    function RHST(arr_size::Int64, irfft_size::Int64)
        new(zeros(arr_size), zeros(irfft_size), zeros(irfft_size))
    end
end



"""
Struct storing the inputs for the function that computes the residuals.
Allows for preallocation of memory.
"""
struct ActionGradFInputsT
    a :: Int64 #rational pair (a, b)
    b :: Int64
    aN :: Int64 #size of fourier expansion in trial function
    nfft :: Int64
    r :: Array{Float64}
    θ :: Array{Float64}
    ζ :: StepRangeLen{Float64}
    prob :: ProblemT
    met :: MetT
    B :: BFieldT
    nlist :: Array{Int64}
    rdot :: RHST #still not really happy with this rhs notation, guess it will have to do.
    θdot :: RHST
    function ActionGradFInputsT(a::Int64, b::Int64, aN::Int64, nfft::Int64, ζ::StepRangeLen{Float64}, prob::ProblemT, met::MetT, B::BFieldT, nlist::UnitRange{Int64})
        #size of the ifft'd arrays.
        arr_size = 2 * aN * nfft
        irfft_size = aN*nfft+1
        rdot = RHST(arr_size, irfft_size)
        θdot = RHST(arr_size, irfft_size)
        new(a, b, aN, nfft, zeros(arr_size), zeros(arr_size), ζ, prob, met, B, nlist, rdot, θdot)
    end
end


"""
Struct storing the derivatives of the 'rhs' data in each equation, includes the computed value and cos/sin coefficients once fft'd.
Used by the jacobian function.
"""
struct DRHST
    drhs :: Array{Float64, 2} 
    dcn :: Array{Float64, 2}
    dsn :: Array{Float64, 2}
    drcn :: Array{Float64, 2}
    drsn :: Array{Float64, 2}
    dθcn :: Array{Float64, 2}
    dθsn :: Array{Float64, 2}
    function DRHST(arr_size::Int64, comb_size::Int64, irfft_size::Int64, aN::Int64)
        new(zeros(arr_size, comb_size), zeros(irfft_size, comb_size), zeros(irfft_size, comb_size), zeros(arr_size, aN+1), zeros(arr_size, aN+1), zeros(arr_size, aN+1), zeros(arr_size, aN+1))
    end
end

"""
Struct storing the inputs for the function that computes the jacobian.
Allows for preallocation of memory.
"""
struct ActionGradJInputsT
    ν_rhs :: Array{Float64} #rhs of equation 3)
    ν_cn :: Array{Float64}
    ν_sn :: Array{Float64}
    cosnza :: Array{Float64, 2}
    sinnza :: Array{Float64, 2}
    drdot :: DRHST
    dθdot :: DRHST
    function ActionGradJInputsT(ip::ActionGradFInputsT)
        arr_size = 2 * ip.aN * ip.nfft
        comb_size = 4*ip.aN+2 #size of array when (most) coefficients are combined for efficent fourier transform.
        irfft_size = ip.aN*ip.nfft+1
        nza = ip.ζ * ip.nlist' / ip.a
        drdot = DRHST(arr_size, comb_size, irfft_size, ip.aN)
        dθdot = DRHST(arr_size, comb_size, irfft_size, ip.aN)
        new(zeros(arr_size), zeros(irfft_size), zeros(irfft_size), cos.(nza), sin.(nza), drdot, dθdot)
    end
end


"""
    action(rational::Tuple{Int64, Int64}, prob::ProblemT, met::MetT, B::BFieldT, M::Int64, N::Int64, sguess=0.5::Float64, nfft=2::Int64)

Extremizes the action for a given rational, returns the fourier coefficients of the solution.
"""
function action(rational::Tuple{Int64, Int64}, prob::ProblemT, met::MetT, B::BFieldT, M::Int64, N::Int64, sguess=0.5::Float64, nfft=2::Int64)


    #a, b is the toroidal, poloidal mode number that defines the rational surfaces a/b
    #surfaces are defined at q=a/b
    #not that b≡p and a≡q in original case, as they use iota.
    a = rational[1]
    b = rational[2]
    #N is the number of fourier points per 2π
    #hence we have 2πa fourier points
    aN = a * N

    #this is a made up additional scaling used by Zhisong, nfft is just an increase to the fourier resolution, which is just scaled arbitrarily by 2.
    MM = 2 * nfft #still very unsure what this is.

    #number of pseudo field lines to be found
    fl = MM * N
    #seems like this should have more to do with M rather than N?
    #fl = N * nfft #number of pseudo field lines to find. Unsure why this has anything to do with N.

    afM = fl * a
    Nfft = afM

    dζ = 2π / fl
    dθ = 2π / afM

    nlist = range(0, aN)
    ζ = range(0, Nfft-1) .* dζ


    #stores the initial guess for all coefficients
    #+2 for n=0 cos terms, +1 for ν equation
    x0 = zeros(4*aN + 2 + 1)

    #these store the minimised values for the coefficeints at each iteration
    rcosarr = zeros(fl , aN+1)
    rsinarr = zeros(fl , aN+1) #does not actually need th +1 (n=0) but kept for symmetry.
    θcosarr = zeros(fl , aN+1)
    θsinarr = zeros(fl , aN+1)
    νarr = zeros(fl)

    #struct storing the inputs needed for action grad
    ipf = ActionGradFInputsT(a, b, aN, nfft, ζ, prob, met, B, nlist)

    #creates struct for storing fft plans and temp arrays
    #for efficient fourier transforming
    ftd = FTDataT(aN, nfft)

    #struct storing the inputs for jacobian function
    ipj = ActionGradJInputsT(ipf)

    #struct storing the coefficients
    coefs = CoefficientsT(aN+1)

    #iterates over different poloidal areas, this creates multiple psuedo fieldlines that minamise the action.
    #these are then combined into a surface below.
    for jpq in 0:fl-1
        #area constraint.
        α = jpq * dθ
        if jpq == 0
            #initial guesses for the alg.
            coefs.rcos[1] = sguess
        else
            #otherwise use the last solution as a starting point.
            coefs.rcos .= rcosarr[jpq, :]
            coefs.rsin .= rsinarr[jpq, :]
            coefs.θcos .= θcosarr[jpq, :]
            coefs.θsin .= θsinarr[jpq, :]
            coefs.ν[1] = νarr[jpq]

            coefs.θcos[1] += dθ
        end
        #combine the values into a single array for the root findning alg's argument
        pack_coeffs!(x0, coefs)

        F0 = similar(x0)

        #creates a standin function that fits the form needed for NLSolve.
        ag_f!(δS, x) = action_grad_f!(δS, x, α, coefs, ipf, ftd)

        #creates a standin function that fits the form needed for NLSolve.
        #this function defines the gradient, improving efficiency.
        ag_j!(JM, x) = action_grad_j!(JM, x, α, coefs, ipf, ftd, ipj)

        #combination of f and j functions, for efficiency
        ag_fj!(δS, JM, x) = action_grad_fj!(δS, JM, x, α, coefs, ipf, ftd, ipj)

        #combined into single object for NLSolve
        df = OnceDifferentiable(ag_f!, ag_j!, ag_fj!, x0, F0)

        #solve for the psuedo field line.
        sol = nlsolve(df, x0)

        #unpacks the single array solution into CoefficientsT struct.
        unpack_coeffs!(sol.zero, coefs, aN+1)

        #stores this field line.
        rcosarr[jpq+1, :] = coefs.rcos
        rsinarr[jpq+1, :] = coefs.rsin
        θcosarr[jpq+1, :] = coefs.θcos
        θsinarr[jpq+1, :] = coefs.θsin
        νarr[jpq+1] = coefs.ν[1]
        
    end

    #wrap_field_lines(rcosarr, rsinarr, θcosarr, θsinarr, MM, M, N, a, b, afM, Nfft)
    #return rcosarr, rsinarr, θcosarr, θsinarr, νarr

    #field lines are wrapped into the qfm surface.
    return wrap_field_lines(rcosarr, rsinarr, θcosarr, θsinarr, MM, M, N, a, b, afM, Nfft)
end



"""

Function that turns the Pseudo field lines into a surface. Probably the least clear function in the entire package.
Names for this are incorrect, will need to be understood one day!
"""
function wrap_field_lines(rcosarr::Array{Float64, 2}, rsinarr::Array{Float64, 2}, θcosarr::Array{Float64, 2}, θsinarr::Array{Float64, 2}, MM::Int64, M::Int64, N::Int64, a::Int64, b::Int64, afM::Int64, Nfft::Int64)
    #inputes are a bit stupid here, probably dont' need all of these.

    #think the params we have here and above would be simpler if we knew what any of them were.
    #this is the number of field lines.
    fM = MM * N 

    #converts the fourier coefficients back to the normal values.
    #this has dimns [number of field lines, whatever the ζarray is, (i think fl*a)!]
    #nfft is hard codes as 2 atm, this whole function needs to be cleaned up.
    r = irfft1D(rcosarr, rsinarr, 2)
    ζ = LinRange(0, 2*a*π, size(r)[end]+1)[1:end-1]

    θcosarr[:, 1] .= 0.0
    θ = irfft1D(θcosarr, θsinarr, 2)


    #so this should be (number of field lines * a, number of coeffs in real space 2*aN)
    #field lines are a bit congfusing atm.
    r2D_alpha = zeros((afM, Nfft))
    θ2D_alpha = zeros((afM, Nfft))
    #display(size(r2D_alpha))

    #this is taken the solutions with α [0, 2π/a] and extending them to [0, 2π]
    #this is done by duplicating the solutions a times,
    #however, the phase of each duplicated point has to be modified.
    #exact way this is done is copied from Zhisong.
    for i in 0:a-1
        idx = mod(b*i, a) #this is the real trick here
        #unsure why this takes the form
        #this dictates the starting point of each 'batch' of field lines
        #i.e. if the periodic orbit has points θ0, θ1, θ2..., θN
        #the idx tells us where the feild lines starting at α0+2π/a should start
        #i.e. they might go θ2, θ3 ..., θN, θ0, θ1.
        #and same for r.
        #this must be forcing the periodicity to be correct, but I cannot replicate in a simpler way.

        r2D_alpha[1+idx * fM : (idx+1)*fM, 1 :(a-i) * N * MM] = r[:, 1+ i * N * MM : end]

        r2D_alpha[1+idx * fM : (idx+1)*fM, 1 + (a-i) * N * MM : end] = r[:, 1: i * N * MM]

        θ2D_alpha[1+idx * fM : (idx+1)*fM, 1 :(a-i) * N * MM] = θ[:, 1+i * N * MM : end]

        θ2D_alpha[1+idx * fM : (idx+1)*fM, 1 + (a-i) * N * MM : end] = θ[:, 1: i * N * MM]
    
    end

    r2D_vartheta = zeros((afM, MM * N))
    θ2D_vartheta = zeros((afM, MM * N))
    #display(size(r2D_vartheta))

    #here we actually get r(ϑ, φ), creating the result we want.
    #this is mapping the evolution in α to evolution in ϑ, zhisongs code kind of explains it
    #but this is also not clear.
    for i in 0:MM * N-1
        #v odd that this is required. but otherwise the mod function removes some of the input???
        arr = 0:afM
        idx = @. mod(arr - i * b, afM) .+ 1

        idx = idx[1:end-1]

        r2D_vartheta[:, i+1] = r2D_alpha[idx, i+1]
        θ2D_vartheta[:, i+1] = θ2D_alpha[idx, i+1]
    end

    #converts the final results back to fourier coefficients, I guess for easier interpolation.
    #this is the only place in the entire method that uses M.
    #very wild, seems like M should be automatically chosen? Given how many conditions there are on it?
    scos_surf, ssin_surf = rfft2D(r2D_vartheta, M, N)
    tcos_surf, tsin_surf = rfft2D(θ2D_vartheta, M, N)

    return scos_surf, tsin_surf, ssin_surf, tcos_surf

end



"""

Function that packs the CoefficientsT struct into a 1d array for solving.
"""
function pack_coeffs!(x::Array{Float64}, coefs::CoefficientsT)

    x .= [coefs.rcos ; coefs.rsin[2:end] ; coefs.θcos ; coefs.θsin[2:end] ; coefs.ν] .+ 1

end


"""

Function that unpacks the 1d list of coefficients into the CoefficientsT struct.
"""
function unpack_coeffs!(x::Array{Float64}, coefs::CoefficientsT, N::Int64)

    coefs.rcos[:] = x[1:N] .- 1
    coefs.rsin[2:end] = x[N+1:2*N-1] .- 1
    coefs.θcos[:] = x[2*N-1+1:3*N-1] .-1
    coefs.θsin[2:end] = x[3*N-1+1:4*N-1-1] .-1
    coefs.ν[1] = x[end] - 1
    coefs.rsin[1] = 0.0
    coefs.θsin[1] = 0.0
    
end
