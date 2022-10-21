#= Code to evolve a quantum-controlled system in time. =#

# include("./utils.jl")
# include("./pulse.jl")
# include("./device.jl")

module Evolutions
using LinearAlgebra
import ..Utils
import ..Pulses
import ..Devices

# ENUMERATIONS TO CONTROL EVOLUTION METHOD
abstract type EvolutionMode end
struct Trotter <: EvolutionMode end
struct Lanczos <: EvolutionMode end



"""
    evolve(
        ψI::AbstractVector{<:Number},
        pulses::AbstractVector{<:Pulses.PulseTemplate},
        device::Devices.Device;
        numsteps::Integer = 2000
    )

Evolve the state `ψI` in time.

The amount of time evolved is determined by the duration of the pulses,
    which are assumed to have equal duration.

# Arguments
- `ψI` initial statevector of `n` qubits each with `nstates` levels
- `pulses` vector of `n` pulse templates
- `device` the `n`-qubit device giving qubit frequencies and couplings
- `numsteps` the number of discrete time units to simulate (ie. Trotter steps)
"""
evolve(
    ψI::AbstractVector{<:Number},
    pulses::AbstractVector{<:Pulses.PulseTemplate},
    device::Devices.Device;
    numsteps::Integer = 2000
) = evolve(ψI, pulses, device, Trotter; numsteps=numsteps)


function evolve(
    ψI::AbstractVector{<:Number},
    # ψI::Vector{Int64},
    pulses::AbstractVector{<:Pulses.PulseTemplate},
    # pulses::Vector{Pulses.BasicSquarePulse},
    device::Devices.Device,
    # device::Devices.Transmon,
    ::Type{Trotter};
    numsteps::Integer = 2000
)
    # INFER NUMBER OF QUBITS AND STATES
    N = length(ψI)                          # SIZE OF STATEVECTOR
    n = length(device)                      # NUMBER OF QUBITS
    nstates = round(Int, N^(1/n))           # NUMBER OF LEVELS ON EACH QUBIT
        # TODO: I feel as though there should be an integer-stable way of doing this...

    # PREPARE TIME DISCRETIZATION
    T = length(pulses[1])                   # TOTAL TIME
        # NOTE: Currently throws error if `pulses` is empty.
        #       We already specify `length(pulses)==length(device)` as a precondition,
        #           but technically `n=0` should maybe be valid..?
        # TODO: Either accommodate in code or forbid in documentation.
    t = range(0,T,numsteps)                 # TIME GRID
        # TODO: My instinct says we skip t=0?
        #       If so, range up to `numsteps+1`
        #           and start from t[2].
    Δt = t[2]-t[1]                          # DURATION OF EACH TROTTER STEP
        # TODO: Control for `numsteps < 2`. Gosh, this section is awkward..!

            # TEMP: Oinam calculates Δt incorrectly, I think. Uncomment below to match his.
            # Δt = T / numsteps

    # CONSTRUCT AND DIAGONALIZE THE DEVICE HAMILTONIAN
    HD = Devices.static_hamiltonian(device, nstates)
    ΛD, UD = eigen(HD)
    UDT = UD'

    # PREPARE CREATION AND ANNIHILATION OPERATORS ON EACH QUBIT, IN DEVICE BASIS
    a1 = Utils.a_matrix(nstates)            # SINGLE-QUBIT ANNIHILATION OPERATOR
    a_ = Vector{Matrix{Number}}(undef, n)   # SINGLE-QUBIT ANNIHILATION FOR EACH QUBIT...
    for q ∈ 1:n                             #   ...but as a multi-qubit operator (`on`)...
        a_[q] = UDT * Utils.on(a1, q, n) * UD   # ...and rotated into the device basis.
    end

    # ROTATE INTO DEVICE BASIS
    ψ = UDT * ψI
    # ψ = I * ψI

    # PERFORM TIME EVOLUTION
    for i ∈ 1:numsteps
        # CONSTRUCT CONTROL HAMILTONIAN (IN DEVICE BASIS)
        HC = zeros(N,N)
        for q ∈ 1:n
            Ω = Pulses.amplitude(pulses[q], t[i])
            ν = Pulses.frequency(pulses[q], t[i])
            T = Ω * exp(im*ν*t[i]) * a_[q]  # RECALL a_[q] IS IN DEVICE BASIS
            HC += (T + T')
        end

        # CONJUGATE WITH ACTION OF (DIAGONALIZED) DEVICE HAMILTONIAN
        exp_itΛD = Diagonal(exp.((im*t[i]) * ΛD))
                # TODO: Any reason to use .* for scalar multiplication? Instinct says no.
        HIC = exp_itΛD * HC * exp_itΛD'     # INTERACTION-PICTURE CONTROL HAMILTONIAN

        # APPLY ACTION OF THE INTERACTION-PICTURE CONTROL HAMILTONIAN
        ψ = exp( (-im*Δt) * HIC) * ψ
                # TODO: Any reason to use .* for scalar multiplication? Instinct says no.
    end

    # ROTATE *OUT* OF DEVICE BASIS
    ψ = UD * ψ
    return ψ

end


function evolve(
    ψI::AbstractVector{<:Number},
    pulses::AbstractVector{<:Pulses.PulseTemplate},
    device::Devices.Device,
    ::Type{Lanczos};
    numsteps::Integer = 2000,
    suzukiorder::Integer = 2,
)
    error("Not yet implemented!")
end

# display(evolve([0,1],[Pulses.BasicSquarePulse(10.0,29.0,[1.0,-1.0],[5])],Devices.Transmon([1],[1],Dict{Devices.QubitCouple,Number}()),Trotter))

end
