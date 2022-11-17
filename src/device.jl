#= Implements encapsulations for device properties, and some standard devices. =#

module Devices

import LinearAlgebra: Hermitian

import ..Utils

# AUXILIARY STRUCTURE TO FACILITATE REPRESENTATION OF QUBIT COUPLINGS
struct QubitCouple
    q1::Int
    q2::Int
    # INNER CONSTRUCTOR: Constrain order so that `Qu...ple(q1,q2) == Qu...ple(q2,q1)`.
    QubitCouple(q1, q2) = q1 > q2 ? new(q2, q1) : new(q1, q2)
end


"""
An abstract representation of a quantum computer describing qubit frequencies and couplings.


# Implementation
A valid Device `device` should implement two methods:
- `Base.length(device)` number of qubits (alternatively, just include an `n` attribute)
- `static_hamiltonian(device, m)` matrix representation with `m` for each qubit
"""
abstract type Device end

"""
    Base.length(device::Device)

Fetch the number of qubits in the device.
"""
Base.length(device::Device) = device.n

"""
    static_hamiltonian(device::Device, m::Int=2)

Construct the matrix representation of the device's static (no pulses) Hamiltonian `H`.

Let ``n``≡`length(device)`:
    `H` acts on the space of ``n`` qubits each with ``m`` levels.
    Thus, the matrix representation is a Hermitian matrix with ``m^n`` rows.

TODO: learn how/when to work with sparse matrices
"""
static_hamiltonian(device::Device, m::Integer=2) = error("Not Implemented")




"""
    Transmon(
        ω::AbstractVector{Float64},
        δ::AbstractVector{Float64},
        gmap::AbstractDict{QubitCouple,Float64}=Dict{QubitCouple,Float64}(),
        n::Int=length(ω),
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
    n::Int
    ω::Vector{Float64}
    δ::Vector{Float64}
    gmap::Dict{QubitCouple,Float64}
    # INNER CONSTRUCTOR: Truncate structures down to n qubits.
    function Transmon(
        ω::Vector{Float64},
        δ::Vector{Float64},
        gmap::Dict{QubitCouple,Float64}=Dict{QubitCouple,Float64}(),
        n::Integer=length(ω),
    )
        # TRUNCATE OFF UNUSED QUBITS FROM FREQUENCY/ANHARMONICITY LISTS
        ω = ω[1:n]
        δ = δ[1:n]
        # FILTER OUT UNUSED QUBITS FROM THE PROVIDING COUPLINGS
        gmap = filter( pair_g -> 1 <= pair_g[1].q1 <= pair_g[1].q2 <= n, gmap)
            # The `filter` closure takes key=>value pairs. `pair_g[1]` is just the key.
        return new(n, ω, δ, gmap)
    end
end

"""
    selectqubits(slice::AbstractVector{Int}, device::Transmon)

Use only a sub-section of a device.

The `slice` vector indicates which qubits to use, and in what order.

"""
function selectqubits(slice::AbstractVector{<:Integer}, device::Transmon)
    # CONSTRUCT PERMUTATION
    σ = zeros(Int, device.n)
    for q in eachindex(slice)
        σ[slice[q]] = q
    end

    # FILTER COUPLING MAP
    gmap = Dict{QubitCouple,Float64}(
        QubitCouple(σ[pair.q1],σ[pair.q2]) => g     # RE-LABEL COUPLINGS BY ORDER IN SLICE
            for (pair, g) in device.gmap
            if all((σ[pair.q1],σ[pair.q2]) .!= 0)   # OMIT COUPLINGS ABSENT FROM SLICE
    )

    return Transmon(device.ω[slice], device.δ[slice], gmap)
end

"""
    static_hamiltonian(device::Transmon, m::Int=2)

Constructs the transmon Hamiltonian
    ``\\sum_q ω_q a_q^† a_q
    - \\sum_q \\frac{δ_q}{2} a_q^† a_q^† a_q a_q
    + \\sum_{⟨pq⟩} g_{pq} (a_p^† a_q + a_q^\\dagger a_p)``.

"""
function static_hamiltonian(device::Transmon, m::Integer=2)
    n = length(device)
    N = m ^ n

    a_ = Utils.algebra(n, m)

    H = zeros(Float64, N,N)

    for q ∈ 1:n
        H += device.ω[q]   * (a_[q]'   * a_[q])     # RESONANCE  TERMS
        H -= device.δ[q]/2 * (a_[q]'^2 * a_[q]^2)   # ANHARMONIC TERMS
    end

    # COUPLING TERMS
    for (pair, g) ∈ device.gmap
        term = g * a_[pair.q1]' * a_[pair.q2]
        H += term + term'
    end

    return Hermitian(H)
end


""" The default four-qubit device used in ctrlq.

I'm not sure exactly where the numbers came from.
The frequencies and anharmonicities are consistent with IBM devices,
    but I don't know of any IBM devices with just four qubits.
Also the anharmonicities and coupling constants form very noticeable patterns...
I'm guessing they were selected somewhat arbitrarily. ^_^
But I'd rather use these arbitrary values than try to come up with my own.

Note however that this device FAILS on ctrlq;
    "dressing" the device eigenbasis DUPLICATES some rows and OMITS others.

"""
DEFAULT_DEVICE = Transmon([
    4.808049015463495,
    4.833254817254613,
    4.940051121317842,
    4.795960998582043,
], [
    0.3101773613134229,
    0.2916170385725456,
    0.3301773613134229,
    0.2616170385725456,
], Dict(
    QubitCouple(1,2) => 0.018312874435769682,
    QubitCouple(1,3) => 0.019312874435769682,
    QubitCouple(1,4) => 0.020312874435769682,
    QubitCouple(2,3) => 0.021312874435769682,
    QubitCouple(2,4) => 0.018312874435769682,
    QubitCouple(3,4) => 0.019312874435769682,
))

end # END MODULE
