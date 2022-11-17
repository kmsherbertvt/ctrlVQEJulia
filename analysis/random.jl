#= Tools to generate random objects for bulk testing. =#


"""
Randomly generate pulses and devices for scalable benchmarking.
"""
module RandomConstructs

import Random: GLOBAL_RNG, AbstractRNG, MersenneTwister
import LinearAlgebra: norm
import Distributions: Gamma, Poisson, Uniform

import ..Devices
import ..Pulses

""" Alternate parameterization of Gamma distribution,
        specifiying mean and standard deviation.
"""
_Gamma(μ, σ) = Gamma((μ/σ)^2, σ^2/μ)

"""
    randomstate([rng=GLOBAL_RNG,] N::Int)

Generate a random quantum state over N complex numbers.

"""
function statevector(rng::AbstractRNG, N::Integer)
    # UNIFORMLY SAMPLE COMPLEX NUMBERS IN QUADRANT 1
    ψ = rand(rng, ComplexF64, N)

    # TRANSFORM TO COVER THE SQUARE INSCRIBED BY UNIT CIRCLE
    ψ .= 2*(ψ .- complex(1,1))

    # NORMALIZE TO UNIT CIRCLE
    ψ ./= norm(ψ)

    #= TODO: These states won't be uniformally distributed.
        Something something "measure theory" "Haar measure" something... =#

    return ψ
end
statevector(N::Integer) = statevector(GLOBAL_RNG, N)



function squarepulse(
    rng::AbstractRNG,   # RANDOM NUMBER GENERATOR
    T::Real,            # MEAN DURATION             (Gamma distribution)
    ν::Real,            # MEAN FREQUENCY            (Gamma distribution)
    Ω::Real,            # MAXIMUM AMPLITUDE (+/-)   (Uniform distribution)
    W::Integer=1;       # MEAN # OF TIME WINDOWS    (Poisson distribution)
    σT=nothing,         # √(σ²) OF T (nothing -> T is fixed)
    σν=nothing,         # √(σ²) OF ν (nothing -> ν is fixed)
    σΩ=nothing,         # nothing -> Ω is fixed, anything else means uniformal sample
    σW=nothing,         # nothing -> W is fixed, anything else means Poisson sample
)
    T = (σT === nothing ? T : rand(rng, _Gamma(T,σT)))
    ν = (σν === nothing ? ν : rand(rng, _Gamma(ν,σν)))
    W = (σW === nothing ? W : rand(rng, Poisson(W)))
    amplitudes = (σΩ === nothing
        ? fill(Ω, W)                    # CONSTANT PULSE
        : rand(rng, Uniform(-Ω,Ω), W)   # BOUNDED PULSE
    )
    steptimes = (σW === nothing
        ? range(0,T,W+1)[2:end-1]       # UNIFORM SPACING
        : sort(T*rand(rng,W-1))         # UNIFORM SAMPLING
    )

    return Pulses.BasicSquarePulse(T, ν, amplitudes, steptimes)
end
squarepulse(
    T::Real, ν::Real, Ω::Real, W::Integer=1; kwargs...
) = squarepulse(GLOBAL_RNG, T, ν, Ω, W; kwargs...)


function transmondevice(
    rng::AbstractRNG,   # RANDOM NUMBER GENERATOR
    n::Integer,         # NUMBER OF QUBITS          (Poisson distribution)
    ω::Real,            # MEAN RESONANT FREQUENCY   (Gamma distribution)
    δ::Real,            # MEAN ANHARMONICITY        (Gamma distribution)
    g::Real = 0.0;      # MEAN COUPLING STRENGTH    (Gamma distribution)
    σC=nothing,         # PAIR-WISE CHANCE OF COUPLING (Bernoulli distribution)
                # nothing -> linear coupling    Vector{QubitCouple} -> specified coupling
    σn=nothing,         # nothing -> n is fixed, anything else means Poisson sample
    σω=ω/9,             # √(σ²) OF RESONANT FREQUENCIES
    σδ=δ/9,             # √(σ²) OF ANHARMONICITY
    σg=g/9,             # √(σ²) OF COUPLING STRENGTH
)
    n = (σn === nothing ? n : rand(rng, Poisson(n)))
    ω = rand(rng, _Gamma(ω, σω), n)
    δ = rand(rng, _Gamma(δ, σδ), n)

    # SELECT COUPLINGS
    couplings = (σC === nothing || σC isa Number
        ? [Devices.QubitCouple(q,q+1) for q in 1:n-1]   # LINEAR COUPLING AT MINIMUM
        : σC                                            # USE GIVEN COUPLINGS
    )
    if σC isa Number                                    # RANDOMLY ADD MORE COUPLINGS
        # ITERATE THROUGH EACH QUBIT PAIR
        for p ∈ 1:n; for q ∈ p+2:n                      # SKIP ADJACENT PAIRS
            # ROLL DIE                                  # IF SUCCESSFUL, ADD QUBIT PAIR
            rand(rng) < σC && push!(couplings, Devices.QubitCouple(p,q))
        end; end
    end
    g = rand(rng, _Gamma(g, σg), length(couplings))     # GENERATE COUPLING STRENGTHS
    gmap = Dict{Devices.QubitCouple,Float64}(couplings .=> g)   # ZIP INTO DICT

    return Devices.Transmon(ω, δ, gmap)
end
transmondevice(
    n::Integer, ω::Real, δ::Real, g::Real = 0.0; kwargs...
) = transmondevice(GLOBAL_RNG, n, ω, δ, g; kwargs...)

end # END MODULE
