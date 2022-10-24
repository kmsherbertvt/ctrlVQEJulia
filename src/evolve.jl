#= Code to evolve a quantum-controlled system in time.

TODO: Just occurred to me, do we need to tell Julia our matrices are Hermitian for eigen()?

=#

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
- `ψI` initial statevector of `n>0` qubits each with `nstates` levels
- `pulses` vector of `n` pulse templates
- `device` the `n`-qubit device giving qubit frequencies and couplings
- `numsteps` the number of discrete time units to simulate (ie. Trotter steps)
             must be positive integer
"""
evolve(
    ψI::AbstractVector{<:Number},
    pulses::AbstractVector{<:Pulses.PulseTemplate},
    device::Devices.Device;
    numsteps::Integer = 2000
) = evolve(ψI, pulses, device, Trotter; numsteps=numsteps)


"""
    evolve(
        ψI::AbstractVector{<:Number},
        pulses::AbstractVector{<:Pulses.PulseTemplate},
        device::Devices.Device,
        ::Type{Trotter};
        numsteps::Integer = 2000
    )

Replicate Trotter method in ctrlq repository,
    which operates in the device basis and treats each step independently.

"""
function evolve(
    ψI::AbstractVector{<:Number},
    pulses::AbstractVector{<:Pulses.PulseTemplate},
    device::Devices.Device,
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
    t = range(0,T,numsteps)                 # TIME GRID
    Δt = numsteps > 1 ? t[2]-t[1] : T       # DURATION OF EACH TROTTER STEP

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
        HIC = exp_itΛD * HC * exp_itΛD'     # INTERACTION-PICTURE CONTROL HAMILTONIAN

        # APPLY ACTION OF THE INTERACTION-PICTURE CONTROL HAMILTONIAN
        ψ = exp( (-im*Δt) * HIC) * ψ
    end

    # ROTATE *OUT* OF DEVICE BASIS
    ψ = UD * ψ
    return ψ

end

"""
    evolve(
        ψI::AbstractVector{<:Number},
        pulses::AbstractVector{<:Pulses.PulseTemplate},
        device::Devices.Device,
        ::Type{Lanczos};
        numsteps::Integer = 2000,
        suzukiorder::Integer = 2,
    )

Apply Lanczos method, combining device action/adjoint into a single repeat step.
    This encourages operating in qubit basis, permitting faster ``H_C`` exponentiation.

"""
function evolve(
    ψI::AbstractVector{<:Number},
    pulses::AbstractVector{<:Pulses.PulseTemplate},
    device::Devices.Device,
    ::Type{Lanczos};
    numsteps::Integer = 2000,
    suzukiorder::Integer = 2,
)
    ######################################################################################
    #                           PRE-SIMULATION SETUP

    # INFER NUMBER OF QUBITS AND STATES
    N = length(ψI)                          # SIZE OF STATEVECTOR
    n = length(device)                      # NUMBER OF QUBITS
    nstates = round(Int, N^(1/n))           # NUMBER OF LEVELS ON EACH QUBIT
        # TODO: I feel as though there should be an integer-stable way of doing this...

    # PREPARE TIME DISCRETIZATION
    T = length(pulses[1])                   # TOTAL TIME
    t = range(0,T,numsteps)                 # TIME GRID
    Δt = numsteps > 1 ? t[2]-t[1] : T       # DURATION OF EACH TROTTER STEP

    # CONSTRUCT AND DIAGONALIZE THE DEVICE HAMILTONIAN
    HD = Devices.static_hamiltonian(device, nstates)
    ΛD, UD = eigen(HD)
    V = UD * Diagonal(exp.((-im*Δt) * ΛD)) * UD'    # REPEATED DEVICE ACTION

    # PREPARE CANONICAL OPERATORS
    a = Utils.a_matrix(nstates)             # SINGLE-QUBIT ANNIHILATION OPERATOR
    Q =      (a + a')                       # CANONICAL COORDINATE OPERATOR
    P = im * (a - a')                       # CANONICAL   MOMENTUM OPERATOR


    # IF SUZUKI MODE IS ACTIVE, PRE-CALCULATE CANONICAL ROTATIONS
    if suzukiorder > 0
        # CONSTRUCT AND DIAGONALIZE CANONICAL OPERATORS
        ΛQq, UQq = eigen(Q)
        ΛPq, UPq = eigen(P)

        # CONSTRUCT PHASE ROTATION MULTIPLIERS (to be scaled by Ω, ν at each time step)
        @assert ΛQq ≈ ΛPq                   # THERE'S ONLY ACTUALLY ONE SET OF EIGENVALUES
        _iΔtΛ = -im * Δt * ΛQq              # THIS WILL DRIVE ALL DYNAMIC TIME EVOLUTIONS

        # CONSTRUCT FULL HILBERT-SPACE CANONICAL ROTATION OPERATORS
        UQ, UP = Matrix(I,1,1), Matrix(I,1,1)
        for q ∈ 1:n
            UQ = kron(UQ, UQq)
            UP = kron(UP, UPq)
        end
    end

    # INTERMEDIATE CONTROL BASIS ROTATIONS (Not needed in lower suzuki orders.)
    UPQ = (suzukiorder >= 1) && UQ' * UP    # ROTATES FROM P -> Q BASIS
    UQP = (suzukiorder >= 2) && UP' * UQ    # ROTATES FROM Q -> P BASIS

    # PREPARE LIGAND OPERATOR TO BIND EACH TIME STEP
    L = (
        (suzukiorder == 0) ?       V        #      APPLY DEVICE ACTION IN QUBIT BASIS
      : (suzukiorder == 1) ? UP' * V * UQ   # Q BASIS -> DEVICE ACTION -> P BASIS
      : (suzukiorder == 2) ? UQ' * V * UQ   # Q BASIS -> DEVICE ACTION -> Q BASIS
      : error("`suzukiorder > 2` is not implemented.")
    )




    ######################################################################################
    #                   SUZUKI-ORDER-SPECIFIC AUXILIARY VARIABLES
    #
    #=
    This part is weird -
        we're using an auxiliary function to summarize the routine for each Trotter step,
        but each Suzuki-order needs a different set of auxiliary variables.
    So, the function takes three vaguely-named variables M1, M2, M3.
        The values they take are set here.

    TODO: a more elegant solution would be to let these auxiliary variables live
        in a global scope. Best is a Dict( nstates => <the thing> ).
        The auxiliary function can just assume that, if it's been called,
            the appropriate dictionaries have been pre-filled.
        Then we can replace `M1, M2, M3` with just `nstates`.
    =#
    M1 = (suzukiorder == 0) ? Q : _iΔtΛ         # Exact solution needs Q and P to build HC.
    M2 = (suzukiorder == 0) ? P : UPQ           # Suzuki needs _iΔtΛ for time evolution,
                                                #   and UPQ to connect the two factors.
    M3 = (suzukiorder == 2) ? UQP : nothing     # suzukiorder=2 connects a third factor...


    ######################################################################################
    #                              BEGIN TIME EVOLUTION

    # FIRST STEP: exp(-𝒊 HD t[0]), but t[0]=0 so this is an identity operation.
    ψ = I * ψI

    # BASIS PRE-ROTATION (Each suzuki order evolves control in a different starting basis.)
    if suzukiorder == 1;    ψ .= UP' * ψ;   end     # PRE-ROTATE INTO P BASIS
    if suzukiorder == 2;    ψ .= UQ' * ψ;   end     # PRE-ROTATE INTO Q BASIS

    # FIRST CONTROL EVOLUTION   (treated separately to give `L` proper "join" behavior)
    _evolvecontrol!(ψ, pulses, t[1], Δt, suzukiorder, M1, M2, M3)

    # TROTTER STEPS
    for i ∈ 2:numsteps
        # CONNECT EACH TIME STEP WITH THE DEVICE ACTION
        ψ .= L * ψ

        # PERFORM CONTROL EVOLUTION
        _evolvecontrol!(ψ, pulses, t[i], Δt, suzukiorder, M1, M2, M3)
    end

    # BASIS POST-ROTATION (Suzuki approximation ends in a different ending basis.)
    if suzukiorder > 0;     ψ .= UQ * ψ;    end     # POST-ROTATE *OUT* OF Q BASIS

    # LAST STEP: exp(𝒊 HD t[numsteps])), ie. exp(-𝒊 HD T)
    ψ .= UD' * ψ                        # ROTATE INTO DEVICE BASIS
    ψ .*= exp.( (im*T) * ΛD)            # ROTATE PHASES FOR ONE LAST TIME EVOLUTION
    ψ .= UD  * ψ                        # ROTATE *OUT* OF DEVICE BASIS

    return ψ
end

"""
    _evolvecontrol!(
        ψ::AbstractVector{<:Number},
        pulses::AbstractVector{<:Pulses.PulseTemplate},
        t, Δt,
        suzukiorder::Integer,
        M1, M2, M3
    )

Auxiliary function to evolve control Hamiltonian in time.

This is a strange function, in that its starting state `ψ`
    and its arguments `M1`, `M2`, `M3` are very different
    depending on value of `suzukiorder`.

`suzukiorder == 0`
    `ψ` is in qubit basis
    `M1` is matrix representation of canonical operator Q ≡    a+a'
    `M2` is    "           "      of     "         "    P ≡ -i(a-a')
    `M3` is unused

`suzukiorder == 1`
    `ψ` is in so-called "P" basis (ie. rotated to diagonal basis of P on each qubit)
    `M1` is eigenvalues of Q operator (or P, they're the same) scaled by -𝒊Δt
    `M2` is unitary matrix to rotate from P -> Q basis
    `M3` is unused

`suzukiorder == 2`
    `ψ` is in so-called "Q" basis (ie. rotated to diagonal basis of Q on each qubit)
    `M1` is eigenvalues of Q operator (or P, they're the same) scaled by -𝒊Δt
    `M2` is unitary matrix to rotate from P -> Q basis
    `M3` is unitary matrix to rotate from Q -> P basis

TODO: The `ψ` will remain weird.
    But `M1`, `M2`, `M3` should be moved to global Dicts, keyed by `nstates`
    (which can be passed as an auxiliary argument, or inferred by `ψ` and `pulses`).

"""
function _evolvecontrol!(
    ψ::AbstractVector{<:Number},
    pulses::AbstractVector{<:Pulses.PulseTemplate},
    t, Δt,
    suzukiorder::Integer,
    M1, M2, M3
)
    n = length(pulses)

    # INTERPRET M1, M2 ACCORDING TO SUZUKI ORDER, AND INITIALIZE FULL-QUBIT-SPACE OPERATORS
    if     suzukiorder == 0
        Q    ,   P, ___ = M1, M2, M3    # CANONICAL OPERATORS
        exp_iΔtHC = Matrix(I,1,1)       # EVOLUTION OF ``H_C``
    else
        _iΔtΛ, UPQ, UQP = M1, M2, M3    # ROTATE FROM Q -> P BASIS
        EQ, EP = ones(1), ones(1)       # EXPONENTIATED PHASES OF Q AND P COMPONENTS
    end

    # CONSTRUCT DIAGONALIZED PHASE ROTATIONS
    #   EXCEPT `suzukiorder==0`, which directly constructs qubit-basis operator.
    for q ∈ 1:n
        # EXTRACT TIME-DEPENDENT COEFFICIENTS
        Ωq = Pulses.amplitude(pulses[q], t)
        νq = Pulses.frequency(pulses[q], t)
        zq = Ωq * exp(im*νq*t)
        xq, yq = real(zq), imag(zq)

        # SUZUKI ORDER 2: EVOLVE Q-COMPONENT ONLY HALF AS FAR (but we'll apply it twice)
        if suzukiorder == 2;    xq /= 2;    end

        # EVOLVE QUBIT IN TIME, AND EXTEND FULL-QUBIT OPERATOR
        if     suzukiorder == 0
            HCq = (xq*Q) + (yq*P)           # SINGLE-QUBIT CONTROL HAMILTONIAN, QUBIT BASIS
            eHCq = exp((-im*Δt) * HCq)      # TIME-EVOLVED
            exp_iΔtHC = kron(exp_iΔtHC, eHCq)   # ATTACHED TO FULL-QUBIT OPERATOR
        else
            EQq = exp.(xq * _iΔtΛ)          # SINGLE-QUBIT TIME-EVOLVED PHASE, Q COMPONENT
            EQ = kron(EQ, EQq )                 # ATTACHED TO FULL-QUBIT OPERATOR

            EPq = exp.(yq * _iΔtΛ)          # SINGLE-QUBIT TIME-EVOLVED PHASE, P COMPONENT
            EP = kron(EP, EPq )                 # ATTACHED TO FULL-QUBIT OPERATOR
        end
    end

    # APPLY ROTATIONS AND PHASES
    if     suzukiorder == 0
        ψ .= exp_iΔtHC * ψ          # APPLY FULL-BODY EVOLUTION OPERATOR
    elseif suzukiorder == 1
        ψ .*= EP                    # ROTATE PHASES FOR TIME EVOLUTION (P BASIS)
        ψ .= UPQ * ψ                # ROTATE FROM P -> Q BASIS
        ψ .*= EQ                    # ROTATE PHASES FOR TIME EVOLUTION (Q BASIS)
    elseif suzukiorder == 2
        ψ .*= EQ                    # ROTATE PHASES FOR TIME EVOLUTION (Q BASIS)
        ψ .= UQP * ψ                # ROTATE FROM Q -> P BASIS
        ψ .*= EP                    # ROTATE PHASES FOR TIME EVOLUTION (P BASIS)
        ψ .= UPQ * ψ                # ROTATE FROM P -> Q BASIS
        ψ .*= EQ                    # ROTATE PHASES FOR TIME EVOLUTION (Q BASIS, 2nd time)
    else
        error("`suzukiorder > 2` is not implemented.")
    end
end

end
