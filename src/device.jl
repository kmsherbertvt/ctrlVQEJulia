#= Implements encapsulations for device properties, and some standard devices. =#

# include("./utils.jl")

module Devices
import ..Utils

# AUXILIARY STRUCTURE TO FACILITATE REPRESENTATION OF QUBIT COUPLINGS
struct QubitCouple
    q1::Integer
    q2::Integer
    # INNER CONSTRUCTOR: Constrain order so that `Qu...ple(q1,q2) == Qu...ple(q2,q1)`.
    QubitCouple(q1, q2) = q1 > q2 ? new(q2, q1) : new(q1, q2)
end


"""
An abstract representation of a quantum computer describing qubit frequencies and couplings.


# Implementation
A valid Device `device` should implement two methods:
- `Base.length(device)` number of qubits (alternatively, just include an `n` attribute)
- `static_hamiltonian(device, nstates)` matrix representation with `nstates` for each qubit
"""
abstract type Device end

"""
    Base.length(device::Device)

Fetch the number of qubits in the device.
"""
Base.length(device::Device) = device.n

"""
    static_hamiltonian(device::Device, nstates::Integer=2)

Construct the matrix representation of the device's static (ie. no fields) Hamiltonian `H`.

Let `n`≡`length(device)`:
    `H` acts on the space of `n` qubits each with `nstates` levels.
    Thus, the matrix representation is a square Matrix with `nstates ^ n` rows.

TODO: learn how/when to work with sparse matrices
"""
static_hamiltonian(device::Device, nstates::Integer=2) = error("Not Implemented")




"""
    Transmon(
        ω::AbstractVector{<:Number},
        δ::AbstractVector{<:Number},
        gmap::AbstractDict{QubitCouple,<:Number},
        n::Integer=length(ω),
    )

A device using transmons as qubits (eg. superconducting quantum computers).

Each transmon is characterized by a resonance frequency between computational states
    and an anharmonicity constant describing higher energy levels.
Device couplings are given as symmetric coupling constants between pairs of qubits.

# Arguments
- `ω` list of `n` resonance frequencies of each qubit.
- `δ` list of `n` anharmonicities of each qubit
- `gmap` keys are qubit pairs (eg. `QubitCouple(1,2)` pairs the qubits indexed 1 and 2),
      and values are static coupling rates
- `n` number of qubits (defaults to length of `ω`)

"""
struct Transmon <: Device
    n::Integer
    ω::AbstractVector{<:Number}
    δ::AbstractVector{<:Number}
    gmap::AbstractDict{QubitCouple,<:Number}
    # INNER CONSTRUCTOR: Truncate structures down to n qubits.
    function Transmon(
        ω::AbstractVector{<:Number},
        δ::AbstractVector{<:Number},
        gmap::AbstractDict{QubitCouple,<:Number},
        n::Integer=length(ω),
    )
        ω = ω[1:n]
        δ = δ[1:n]
        gmap = Dict(pair=>g for (pair, g) ∈ gmap if 1 <= pair.q1 <= pair.q2 <= n)
        return new(n, ω, δ, gmap)
    end
end
# TODO: Mayhaps the inner constructor should be strictly four arguments,
#       and we make an outer constructor that infers n?

"""
    static_hamiltonian(device::Transmon, nstates::Integer)

Constructs the transmon Hamiltonian
    ``\\sum_q ω_q a_q^† a_q
    - \\sum_q \\frac{δ_q}{2} a_q^† a_q^† a_q a_q
    + \\sum_{⟨pq⟩} g_{pq} (a_p^† a_q + a_q^\\dagger a_p)``.

NOTE: This function is not implemented very efficiently.
    But, it only needs to happen once...
"""
function static_hamiltonian(device::Transmon, nstates::Integer=2)
    n = length(device)
    N = nstates ^ n

    a_ = Utils.a_matrix(nstates)
    aT = a_'

    H = zeros(N,N)

    for q ∈ 1:n
        a_q = Utils.on(a_,q,n)
        aTq = Utils.on(aT,q,n)
        H += device.ω[q]   * (aTq * a_q)        # RESONANCE  TERMS
        H -= device.δ[q]/2 * (aTq^2 * a_q^2)    # ANHARMONIC TERMS
    end

    # COUPLING TERMS
    for (pair, g) ∈ device.gmap
        a_1 = Utils.on(a_,pair.q1,n)
        a_2 = Utils.on(a_,pair.q2,n)
        aT1 = Utils.on(aT,pair.q1,n)
        aT2 = Utils.on(aT,pair.q2,n)

        H += g * (aT1 * a_2 + aT2 * a_1)
    end

    return H
end









# CONSTRUCT A DEFAULT_DEVICE FOR CONVENIENT TESTING
DEFAULT_DEVICE = Transmon(
    2π .* [4.8080490154634950, 4.8332548172546130],
    2π .* [0.3101773613134229, 0.2916170385725456],
    Dict(
        QubitCouple(1,2) => 2π*0.018312874435769682,
    ),
)

end
