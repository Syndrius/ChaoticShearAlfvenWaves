#this may not be the best use of time.
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
struct ActionT
    a :: Int64
    b :: Int64
    qN :: Int64
    nfft :: Int64
    r :: Array{Float64}
    θ :: Array{Float64}
    ζ :: StepRangeLen{Float64}
    prob :: ProblemT
    met :: MetT
    B :: BFieldT
    nlist :: Array{Int64}
    rrhs :: Array{Float64}
    θrhs :: Array{Float64}
    function ActionT(a::Int64, b::Int64, qN::Int64, nfft::Int64, ζ::StepRangeLen{Float64}, prob::ProblemT, met::MetT, B::BFieldT, nlist::UnitRange{Int64})
        #perhaps the inputs are a bit restrictive...
        #fft_size = qN * nfft + 1
        ifft_size = 2 * qN * nfft
        new(a, b, qN, nfft, zeros(ifft_size), zeros(ifft_size), ζ, prob, met, B, nlist, zeros(ifft_size), zeros(ifft_size))
    end
end

function action(a::Int64, b::Int64, prob::ProblemT, met::MetT, B::BfieldT, M::Int64, N::Int64, sguess=0.5::Float64, nfft=2::Int64)
    #a, b is the toroidal, poloidal mode number that defines the rational surfaces a/b
    #surfies is defined at q=a/b
    #not that b≡p and a≡q in original case, as they use iota.

    #N is the number of fourier harmonics included in the trial function.
    #this gets expanded to q*N, as the new field line will have 2πq periodicity
    #so the domain is expanded, requiring more fourier harmonics.
    qN = q * N

    MM = 2 * nfft #still very unsure what this is.

    #number of pseudo field lines to be found
    fl = MM * N

    #r(ζ) = r_0^c + ∑_{n=1}^qN (r_n^c cos(nζ/q) + r_n^s sin(nζ/q))
    #θ(ζ) = θ_0^c + ι*ζ + ∑_{n=1}^qN (θ_n^c cos(nζ/q) + θ_n^s sin(nζ/q))

    nlist = range(0, qN)

    #array storing the unknowns
    #[[ν], [rcos], [θsin[2:end]], [rsin[2:end]], [θcos]]
    #sin[1] components are ignored as sin(0) = 0.
    #qN counts from n=1 to n=qN
    #4 * qN is the 4 set of coeffs for n=1 to n=qN
    #+1 for ν
    #+2 for the two cos terms which are n=0

    x = zeros(4*qN + 1 + 2)

    #these store the minimised values for the coefficeints at each iteration
    rcosarr = zeros(fl , qN+1)
    rsinarr = zeros(fl , qN+1) #does not actually need th +1 (n=0) but kept for symmetry.
    θcosarr = zeros(fl , qN+1)
    θsinarr = zeros(fl , qN+1)
    nvarr = zeros(fl)

    #size of arrays when fourier transformed, scaled by nfft to prevent aliasing
    #+1 is for n=0 terms
    ip = ActionT(a, b, qN, fft_size, 


#this one computes the action gradient and the second derivative.
#this will also hopefully be a bit clearer.
function comb_action_grad!()
   
    unpack_coeffs!(x, coefs, ip.Ntor)

    #should be changed to b / a to match other work.
    iota = ip.p / ip.q

    #note this is the area of the trial function
    area = coefs.θcos[1]


    #inverse fourier transform, done for r and θ at the same time.
    get_r_t!(ip.r, ip.θ, coefs, ftd.ift_1D_p, ftd.ift_r1D, ftd.ift_θ1D, ip.Ntor)

    #reflects that there is an additional term in the trial function for θ
    ip.θ .+= iota .* ip.ζ

    for i in 1:1:length(ip.r)
        ip.prob.met(ip.met, ip.r[i], ip.θ[i], ip.ζ[i], ip.prob.geo.R0)

        compute_B!(ip.B, ip.met, ip.prob.q, ip.prob.isls, ip.r[i], ip.θ[i], ip.ζ[i])

        #these are awful names,
        #θdot is what θdot is equal to, when the gradient is zero.
        #this is why Zhisong used rhs terminology.
        ip.θdot[i] = ip.B.B[2] / ip.B.B[3]

        ip.rdot[i] = ip.B.B[1] / ip.B.B[3] - coefs.nv[1] / (ip.B.B[3] * ip.met.J[1])


