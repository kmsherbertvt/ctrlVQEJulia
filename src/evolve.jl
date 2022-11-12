#= Code to evolve a quantum-controlled system in time. =#

module Evolutions

import LinearAlgebra: eigen, Hermitian, Diagonal
import DifferentialEquations: ODEProblem, solve
import KrylovKit: exponentiate
import TensorOperations: ncon

import ..Utils
import ..Pulses
import ..Devices

# ENUMERATIONS TO CONTROL EVOLUTION METHOD
abstract type EvolutionMode end
struct ODE     <: EvolutionMode end
struct Direct  <: EvolutionMode end
struct Lanczos <: EvolutionMode end
struct Rotate  <: EvolutionMode end
struct Prediag <: EvolutionMode end

# ENUMERATIONS TO CONTROL INPUT/OUTPUT BASIS
abstract type IOBasisMode end
struct QubitBasis  <: IOBasisMode end
struct DeviceBasis <: IOBasisMode end

# ENUMERATIONS TO CONTROL QUBIT-WISE OPERATIONS
abstract type QubitApplyMode end
struct Kronec <: QubitApplyMode end
struct Tensor <: QubitApplyMode end



"""
    evolve(ψI, args...; kwargs...)

Shorthand for a not-in-place function to evolve a state `ψI` in time.

This just copies ψI to a new variable then calls `evolve!`, so find documentation there!

"""
function evolve(ψI, args...; kwargs...)
    ψ = copy(ψI)
    evolve!(ψ, args...; kwargs...)
    return ψ
end


"""
    evolve!(
        ψ::AbstractVector{<:Number},
        pulses::AbstractVector{<:Pulses.PulseTemplate},
        device::Devices.Device;
        iobasis::IOBasisMode = DeviceBasis()
    )

Evolve the state `ψ` in time.

# Arguments
- `ψ` initial statevector of `n>0` qubits each with `m` levels
- `pulses` vector of `n` pulse templates
- `device` the `n`-qubit device giving qubit frequencies and couplings
- `iobasis` which basis the input is interpreted as, and output converted to

The amount of time evolved is determined by the duration of the pulses,
    which are assumed to have equal duration.

All evolutions are performed in the interaction-picture,
    meaning the static time dependence due to device eigenenergies is implicit.
In other words, if there are no pulses, the input state doesn't change at all.
This is distinct from the textbook "Schrodinger picture",
    where an eigenvector of the device Hamiltonian with eigenvalue ``ω``
    would after a time ``T`` incur a phase shift ``\\exp(𝒊ωT)``.

The input/output basis may be controlled with the `iobasis` keyword argment.
The choices are:
- `QubitBasis()` computational states correspond to product states of individual qubits.
- `DeviceBasis()` computational states are eigenstates of the static device Hamiltonian.
  Note that the default device Hamiltonian factorization uses a so-called  "dressed basis"
    (see `Utils.dressedbasis`) but alternative factorizations can usually be passed in
    manually as keyword arguments when calling `evolve!` with a specific algorithm.

The specific algorithm can be selected by passsing an `EvolutionMode`
    as a final positional argument.
These are enumeration structs defined within the `Evolutions` module;
    study the method headers below for supported types.
Note you must pass the raw type (eg. `ODE`) rather than a singleton object (eg. `ODE()`).

Each algorithm includes additional keyword arguments;
    consult the individual method documentation for details.
In particular, note that most algorithms involve a "setup" phase,
    doing expensive calculations that only need to be done once ever for a given device.
So, if you are planning to call `evolve!` *more* than once (eg. optimizing pulse parameters)
    you should take advantage of the optional keyword arguments
    that let you skip those expensive calculations.
Study the method headers for each mode in the code for the requisite details.

"""
evolve!(
    ψ::AbstractVector{<:Number},
    pulses::AbstractVector{<:Pulses.PulseTemplate},
    device::Devices.Device;
    iobasis::IOBasisMode = DeviceBasis()
) = evolve!(ψ, pulses, device, Rotate; iobasis=iobasis)




"""
    evolve!(
        ψ::AbstractVector{<:Number},
        pulses::AbstractVector{<:Pulses.PulseTemplate},
        device::Devices.Device,
        ::Type{ODE};
        iobasis::IOBasisMode = DeviceBasis(),

        # INFERRED VALUES (relatively fast, but pass them in to minimize allocations)
        N = length(ψ),                  # SIZE OF STATEVECTOR
        n = length(device),             # NUMBER OF QUBITS
        m = round(Int, N^(1/n)),        # NUMBER OF LEVELS ON EACH QUBIT
        T = length(pulses[1]),          # TOTAL DURATION OF EVOLUTION

        # CALCULATED VALUES (pass these in to significantly speed up optimizations)
        ΛD = nothing,                   # EIGENVALUES OF STATIC HAMILTONIAN
        UD = nothing,                   # CORRESPONDING EIGENVECTORS
        a_ = nothing,                   # LIST OF ANNIHILATION OPERATORS, IN STATIC BASIS
    )

Numerically integrate Schrodinger's equation.

This method uses state-of-the-art algorithms to adaptively integrate
    a set of coupled differential equations.
That is my way of saying I don't actually know how this method works.

It is extremely accurate for a single qubit,
    so I consider it the "gold standard" for accuracy estimations.
I've no idea how its runtime scales with system size.

"""
function evolve!(
    ψ::AbstractVector{<:Number},
    pulses::AbstractVector{<:Pulses.PulseTemplate},
    device::Devices.Device,
    ::Type{ODE};
    iobasis::IOBasisMode = DeviceBasis(),

    # INFERRED VALUES (relatively fast, but pass them in to minimize allocations)
    N = length(ψ),                  # SIZE OF STATEVECTOR
    n = length(device),             # NUMBER OF QUBITS
    m = round(Int, N^(1/n)),        # NUMBER OF LEVELS ON EACH QUBIT
    T = length(pulses[1]),          # TOTAL DURATION OF EVOLUTION

    # CALCULATED VALUES (pass these in to significantly speed up optimizations)
    ΛD = nothing,                   # EIGENVALUES OF STATIC HAMILTONIAN
    UD = nothing,                   # CORRESPONDING EIGENVECTORS
    a_ = nothing,                   # LIST OF ANNIHILATION OPERATORS, IN STATIC BASIS
)
    ######################################################################################
    #                            PRELIMINARY CALCULATIONS

    if any((ΛD, UD) .=== nothing)
        HD = Devices.static_hamiltonian(device, m)  # DEVICE HAMILTONIAN
        ΛD, UD = Utils.dressedbasis(HD)             # DEVICE EIGENVALUES AND EIGENVECTORS
    end
    if a_ === nothing
        a_ = Utils.algebra(n, m, basis=UD)          # LIST OF ROTATED ANNIHILATION OPERATORS
    end

    if iobasis isa QubitBasis;  ψ .= UD' * ψ;   end;    # ROTATE INTO DEVICE BASIS

    ######################################################################################
    #                       DEFINE AND SOLVE DIFFERENTIAL EQUATIONS

    # DEFINE THE INTERACTION-PICTURE HAMILTONIAN FOR A GIVEN TIME
    function interaction!(du, u, p, t)
        # CONSTRUCT CONTROL HAMILTONIAN (IN DEVICE BASIS)
        HC = zeros(N,N)
        for q ∈ 1:n
            Ω = Pulses.amplitude(pulses[q], t)
            ν = Pulses.frequency(pulses[q], t)
            z = Ω * exp(im*ν*t)
            HC += z*a_[q] + z'*a_[q]'
        end

        # CONJUGATE WITH ACTION OF (DIAGONALIZED) DEVICE HAMILTONIAN
        expD = Diagonal(exp.((im*t) * ΛD))  # DEVICE ACTION
        HIC = expD * HC * expD'     # INTERACTION-PICTURE CONTROL HAMILTONIAN

        # TODO: pre-allocate HC, expD, and HIC into p, perhaps?

        du .= -im * HIC * u
    end

    # SOLVE THE SYSTEM
    #= NOTE:
        This method autoselects an algorithm based on I have no idea what,
            meaning I have no idea what the time complexity or accuracy are likely to be,
            or how I should expect them to scale with increasing system size.
        But, it works *pretty* well for the single-qubit case,
            so I'm willing to treat it as a sort of black-box standard.
    =#
    schrodinger = ODEProblem(interaction!, ψ, (0.0, T))
    solution = solve(schrodinger, save_everystep=false)     # TIME-CONSUMING STEP

    # WRITE FINAL SOLUTION TO THE GIVEN STATEVECTOR
    ψ .= solution.u[end]

    # RE-NORMALIZE THIS STATE
    ψ .= ψ / √abs(ψ'*ψ)

    ######################################################################################

    if iobasis isa QubitBasis;  ψ .= UD  * ψ;   end;    # ROTATE *OUT* OF DEVICE BASIS
end





"""
    evolve!(
        ψ::AbstractVector{<:Number},
        pulses::AbstractVector{<:Pulses.PulseTemplate},
        device::Devices.Device,
        ::Type{Direct};
        iobasis::IOBasisMode = DeviceBasis(),
        numsteps::Integer = 2000,

        # INFERRED VALUES (relatively fast, but pass them in to minimize allocations)
        N = length(ψ),                  # SIZE OF STATEVECTOR
        n = length(device),             # NUMBER OF QUBITS
        m = round(Int, N^(1/n)),        # NUMBER OF LEVELS ON EACH QUBIT
        T = length(pulses[1]),          # TOTAL DURATION OF EVOLUTION
        t_= range(0,T,numsteps+1),      # TIME GRID
        Δt= T / numsteps,               # DURATION OF EACH TIME STEP

        # CALCULATED VALUES (pass these in to significantly speed up optimizations)
        ΛD = nothing,                   # EIGENVALUES OF STATIC HAMILTONIAN
        UD = nothing,                   # CORRESPONDING EIGENVECTORS
        a_ = nothing,                   # LIST OF ANNIHILATION OPERATORS, IN STATIC BASIS
    )

Trotterize the time-evolution operator,
    directly exponentiating ``\\exp(-𝒊·Δt·H)`` at each time step.

The keyword argument `numsteps` specifies the number of time steps;
    time scales linearly, and accuracy scales inversely.

This method is the only one that invokes an ``O(N^3)`` operation at every time step.
Don't use it, except to illustrate how much better other methods are. ^_^

"""
function evolve!(
    ψ::AbstractVector{<:Number},
    pulses::AbstractVector{<:Pulses.PulseTemplate},
    device::Devices.Device,
    ::Type{Direct};
    iobasis::IOBasisMode = DeviceBasis(),
    numsteps::Integer = 2000,

    # INFERRED VALUES (relatively fast, but pass them in to minimize allocations)
    N = length(ψ),                  # SIZE OF STATEVECTOR
    n = length(device),             # NUMBER OF QUBITS
    m = round(Int, N^(1/n)),        # NUMBER OF LEVELS ON EACH QUBIT
    T = length(pulses[1]),          # TOTAL DURATION OF EVOLUTION
    t_= range(0,T,numsteps+1),      # TIME GRID
    Δt= T / numsteps,               # DURATION OF EACH TIME STEP

    # CALCULATED VALUES (pass these in to significantly speed up optimizations)
    ΛD = nothing,                   # EIGENVALUES OF STATIC HAMILTONIAN
    UD = nothing,                   # CORRESPONDING EIGENVECTORS
    a_ = nothing,                   # LIST OF ANNIHILATION OPERATORS, IN STATIC BASIS
)
    ######################################################################################
    #                            PRELIMINARY CALCULATIONS

    if any((ΛD, UD) .=== nothing)
        HD = Devices.static_hamiltonian(device, m)  # DEVICE HAMILTONIAN
        ΛD, UD = Utils.dressedbasis(HD)             # DEVICE EIGENVALUES AND EIGENVECTORS
    end; if a_ === nothing
        a_ = Utils.algebra(n, m, basis=UD)          # LIST OF ROTATED ANNIHILATION OPERATORS
    end

    if iobasis isa QubitBasis;  ψ .= UD' * ψ;   end;    # ROTATE INTO DEVICE BASIS

    ######################################################################################
    #                                 TIME EVOLUTION

    # FIRST TIME STEP   (use Δt/2 for first and last time step)
    ψ .= exp( (-im*Δt/2) * _interactionhamiltonian(pulses, ΛD, a_, t_[1]; N=N, n=n)) * ψ

    for i ∈ 2:numsteps
        ψ .= exp( (-im*Δt) * _interactionhamiltonian(pulses, ΛD, a_, t_[i]; N=N, n=n)) * ψ
    end

    # LAST TIME STEP    (use Δt/2 for first and last time step)
    ψ .= exp( (-im*Δt/2) * _interactionhamiltonian(pulses, ΛD, a_, t_[end]; N=N, n=n)) * ψ

    ######################################################################################

    # RE-NORMALIZE THIS STATE
    ψ .= ψ / √abs(ψ'*ψ)

    if iobasis isa QubitBasis;  ψ .= UD  * ψ;   end;    # ROTATE *OUT* OF DEVICE BASIS
end


"""
    evolve!(
        ψ::AbstractVector{<:Number},
        pulses::AbstractVector{<:Pulses.PulseTemplate},
        device::Devices.Device,
        ::Type{Direct};
        iobasis::IOBasisMode = DeviceBasis(),
        numsteps::Integer = 2000,

        # INFERRED VALUES (relatively fast, but pass them in to minimize allocations)
        N = length(ψ),                  # SIZE OF STATEVECTOR
        n = length(device),             # NUMBER OF QUBITS
        m = round(Int, N^(1/n)),        # NUMBER OF LEVELS ON EACH QUBIT
        T = length(pulses[1]),          # TOTAL DURATION OF EVOLUTION
        t_= range(0,T,numsteps+1),      # TIME GRID
        Δt= T / numsteps,               # DURATION OF EACH TIME STEP

        # CALCULATED VALUES (pass these in to significantly speed up optimizations)
        ΛD = nothing,                   # EIGENVALUES OF STATIC HAMILTONIAN
        UD = nothing,                   # CORRESPONDING EIGENVECTORS
        a_ = nothing,                   # LIST OF ANNIHILATION OPERATORS, IN STATIC BASIS
    )

Trotterize the time-evolution operator,
    calculating the matrix exponential action ``\\exp(-𝒊·Δt·H) |ψ⟩`` at each time step.

The keyword argument `numsteps` specifies the number of time steps;
    time scales linearly, and accuracy scales inversely.

This method is the only one that invokes an ``O(N^3)`` operation at every time step.
Don't use it, except to illustrate how much better other methods are. ^_^

"""
function evolve!(
    ψ::AbstractVector{<:Number},
    pulses::AbstractVector{<:Pulses.PulseTemplate},
    device::Devices.Device,
    ::Type{Lanczos};
    iobasis::IOBasisMode = DeviceBasis(),
    numsteps::Integer = 2000,

    # INFERRED VALUES (relatively fast, but pass them in to minimize allocations)
    N = length(ψ),                  # SIZE OF STATEVECTOR
    n = length(device),             # NUMBER OF QUBITS
    m = round(Int, N^(1/n)),        # NUMBER OF LEVELS ON EACH QUBIT
    T = length(pulses[1]),          # TOTAL DURATION OF EVOLUTION
    t_= range(0,T,numsteps+1),      # TIME GRID
    Δt= T / numsteps,               # DURATION OF EACH TIME STEP

    # CALCULATED VALUES (pass these in to significantly speed up optimizations)
    ΛD = nothing,                   # EIGENVALUES OF STATIC HAMILTONIAN
    UD = nothing,                   # CORRESPONDING EIGENVECTORS
    a_ = nothing,                   # LIST OF ANNIHILATION OPERATORS, IN STATIC BASIS
)
    ######################################################################################
    #                            PRELIMINARY CALCULATIONS

    if any((ΛD, UD) .=== nothing)
        HD = Devices.static_hamiltonian(device, m)  # DEVICE HAMILTONIAN
        ΛD, UD = Utils.dressedbasis(HD)             # DEVICE EIGENVALUES AND EIGENVECTORS
    end; if a_ === nothing
        a_ = Utils.algebra(n, m, basis=UD)          # LIST OF ROTATED ANNIHILATION OPERATORS
    end

    if iobasis isa QubitBasis;  ψ .= UD' * ψ;   end;    # ROTATE INTO DEVICE BASIS

    ######################################################################################
    #                                 TIME EVOLUTION

    # FIRST TIME STEP   (use Δt/2 for first and last time step)
    ψ .= exponentiate(
        _interactionhamiltonian(pulses, ΛD, a_, t_[1]; N=N, n=n), -im*Δt/2,  ψ
    )[1]        # `exponentiate` RETURNS A TUPLE, WE CARE ONLY ABOUT FIRST ELEMENT

    for i ∈ 2:numsteps
        ψ .= exponentiate(
            _interactionhamiltonian(pulses, ΛD, a_, t_[i]; N=N, n=n), -im*Δt, ψ
        )[1]    # `exponentiate` RETURNS A TUPLE, WE CARE ONLY ABOUT FIRST ELEMENT
    end

    # LAST TIME STEP    (use Δt/2 for first and last time step)
    ψ .= exponentiate(
        _interactionhamiltonian(pulses, ΛD, a_, t_[end]; N=N, n=n), -im*Δt/2, ψ
    )[1]        # `exponentiate` RETURNS A TUPLE, WE CARE ONLY ABOUT FIRST ELEMENT

    ######################################################################################

    # RE-NORMALIZE THIS STATE
    ψ .= ψ / √abs(ψ'*ψ)

    if iobasis isa QubitBasis;  ψ .= UD  * ψ;   end;    # ROTATE *OUT* OF DEVICE BASIS
end

"""
    _interactionhamiltonian(
        pulses::AbstractVector{<:Pulses.PulseTemplate}, # PULSE TEMPLATES FOR EACH QUBIT
        ΛD::AbstractVector{<:Number},                   # EIGENVALUES OF STATIC HAMILTONIAN
        a_::AbstractVector{<:AbstractMatrix{<:Number}}  # LIST OF ROTATED ANNIHILATION OPS
        t::Number;                                      # TIME POINT

        # INFERRED VALUES (relatively fast, but pass them in to minimize allocations)
        N = length(ΛD),                                 # SIZE OF STATEVECTOR
        n = length(pulses),                             # NUMBER OF QUBITS
    )

Construct the interaction-picture Hamiltonian for a given time point.

Mathematically, this is ``\\exp(𝒊·t·H)·V(t)·exp(-𝒊·t·H)``,
    where H is static device Hamiltonian and
    where V(t) is the control Hamiltonian
        ``\\sum_q Ω_q(t)[\\exp(𝒊ν_qt) a_q + \\exp(-𝒊ν_qt) a_q^\\dagger]``
Computationally, assume we're in the device basis so the conjugating factor is diagonal.

"""
function _interactionhamiltonian(
    pulses::AbstractVector{<:Pulses.PulseTemplate}, # PULSE TEMPLATES FOR EACH QUBIT
    ΛD::AbstractVector{<:Number},                   # NUMBER OF LEVELS ON EACH QUBIT
    a_::AbstractVector{<:AbstractMatrix{<:Number}}, # LIST OF ROTATED ANNIHILATION OPS
    t::Number;                                      # TIME POINT

    # INFERRED VALUES (relatively fast, but pass them in to minimize allocations)
    N = length(ΛD),                                 # SIZE OF STATEVECTOR
    n = length(pulses),                             # NUMBER OF QUBITS
)
    # CONSTRUCT CONTROL HAMILTONIAN (IN DEVICE BASIS)
    HC = zeros(N,N)
    for q ∈ 1:n
        Ω = Pulses.amplitude(pulses[q], t)
        ν = Pulses.frequency(pulses[q], t)
        z = Ω * exp(im*ν*t)
        HC += z*a_[q] + z'*a_[q]'
    end

    # CONJUGATE WITH ACTION OF (DIAGONALIZED) DEVICE HAMILTONIAN
    expD = Diagonal(exp.((im*t) * ΛD))  # DEVICE ACTION
    HIC = expD * HC * expD'     # INTERACTION-PICTURE CONTROL HAMILTONIAN

    return Hermitian(HIC)
    # TODO: pre-allocate HC and expD.
end





"""
    evolve!(
        ψ::AbstractVector{<:Number},
        pulses::AbstractVector{<:Pulses.PulseTemplate},
        device::Devices.Device,
        ::Type{Rotate};
        iobasis::IOBasisMode = DeviceBasis(),
        numsteps::Integer = 2000,
        qubitapplymode::QubitApplyMode = Kronec(),

        # INFERRED VALUES (relatively fast, but pass them in to minimize allocations)
        N = length(ψ),                      # SIZE OF STATEVECTOR
        n = length(device),                 # NUMBER OF QUBITS
        m = round(Int, N^(1/n)),            # NUMBER OF LEVELS ON EACH QUBIT
        T = length(pulses[1]),              # TOTAL DURATION OF EVOLUTION
        t_= range(0,T,numsteps+1),          # TIME GRID
        Δt= T / numsteps,                   # DURATION OF EACH TIME STEP

        # CALCULATED VALUES (pass these in to significantly speed up optimizations)
        ΛD = nothing,                       # EIGENVALUES OF STATIC HAMILTONIAN
        UD = nothing,                       # CORRESPONDING EIGENVECTORS
        V  = nothing,                       # REPEATED DEVICE ACTION
        a  = nothing,                       # SINGLE-QUBIT ANNIHILATION OPERATOR

        # PRE-ALLOCATIONS (for those that want every last drop of efficiency...)
        O_ = [Matrix{ComplexF64}(undef, m,m) for q ∈ 1:n]   # QUBIT-WISE DRIVE OPERATORS
    )

Trotterize the time evolution operator,
    but switch back and forth between the static and drive bases
    so that the time evolution in each is more efficient.

This method invokes a number of matrix-vector multiplications (``O(N^2)``),
    and some small matrix exponentiations (``O(n·m^3)``) at each time step.
Of the more efficient algorithms, this is the easiest to understand,
    and it will tend to be the most performant also.

The keyword argument `numsteps` specifies the number of time steps;
    time scales linearly, and accuracy scales inversely.

The keyword argument `qubitapplymode` specifies the numerical approach to
    applying the qubit-wise drive operators at each time step.
The choices are:
- `Kronec()` this method combines all qubit-wise operators into a single N×N matrix,
    then applies them with a single matrix-vector multiplication.
- `Tensor()` this method reshapes the statevector into an n-dimensional array,
    and performs a tensor contraction over each qubit-wise operator.
  In principle this one should scale significantly better than `Kronec`,
    but in practice the overhead from tensor manipulation may be steep.
TODO: The horrible overhead on tensor algebra might just vanish by caching...

"""
function evolve!(
    ψ::AbstractVector{<:Number},
    pulses::AbstractVector{<:Pulses.PulseTemplate},
    device::Devices.Device,
    ::Type{Rotate};
    iobasis::IOBasisMode = DeviceBasis(),
    numsteps::Integer = 2000,
    qubitapplymode::QubitApplyMode = Kronec(),

    # INFERRED VALUES (relatively fast, but pass them in to minimize allocations)
    N = length(ψ),                      # SIZE OF STATEVECTOR
    n = length(device),                 # NUMBER OF QUBITS
    m = round(Int, N^(1/n)),            # NUMBER OF LEVELS ON EACH QUBIT
    T = length(pulses[1]),              # TOTAL DURATION OF EVOLUTION
    t_= range(0,T,numsteps+1),          # TIME GRID
    Δt= T / numsteps,                   # DURATION OF EACH TIME STEP

    # CALCULATED VALUES (pass these in to significantly speed up optimizations)
    ΛD = nothing,                       # EIGENVALUES OF STATIC HAMILTONIAN
    UD = nothing,                       # CORRESPONDING EIGENVECTORS
    V  = nothing,                       # REPEATED DEVICE ACTION
    a  = nothing,                       # SINGLE-QUBIT ANNIHILATION OPERATOR

    # PRE-ALLOCATIONS (for those that want every last drop of efficiency...)
    O_ = [Matrix{ComplexF64}(undef, m,m) for q ∈ 1:n]   # QUBIT-WISE DRIVE OPERATORS
)
    ######################################################################################
    #                            PRELIMINARY CALCULATIONS

    if any((ΛD, UD) .=== nothing)
        HD = Devices.static_hamiltonian(device, m)  # DEVICE HAMILTONIAN
        ΛD, UD = Utils.dressedbasis(HD)             # DEVICE EIGENVALUES AND EIGENVECTORS
    end; if V === nothing
        V = UD* Diagonal(exp.((-im*Δt) * ΛD)) *UD'  # REPEATED DEVICE ACTION
    end; if a === nothing
        a = Utils.a_matrix(m)                       # SINGLE-QUBIT ANNIHILATION OPERATOR
    end

    if iobasis isa DeviceBasis; ψ .= UD * ψ;    end;    # ROTATE *OUT* OF DEVICE BASIS

    ######################################################################################
    #                                 TIME EVOLUTION

    #= NOTE: The very first step is, mathematically, exp(-𝒊 HD t_[1]),
        but since t_[1]=0, this is an identity operation and we can skip it. =#

    # APPLY FIRST PULSE DRIVES  (use Δt/2 for first and last time step)
    _preparequbitdrives(pulses, m, t_[1], Δt/2; n=n, a=a, O_=O_)
    _applyqubitoperators!(ψ, O_, qubitapplymode; N=N, n=n, m=m)

    for i ∈ 2:numsteps
        # CONNECT EACH TIME STEP WITH THE DEVICE ACTION
        ψ .= V * ψ

        # APPLY PULSE DRIVES
        _preparequbitdrives(pulses, m, t_[i], Δt; n=n, a=a, O_=O_)
        _applyqubitoperators!(ψ, O_, qubitapplymode; N=N, n=n, m=m)
    end

    # APPLY LAST PULSE DRIVES   (use Δt/2 for first and last time step)
    ψ .= V * ψ
    _preparequbitdrives(pulses, m, t_[end], Δt/2; n=n, a=a, O_=O_)
    _applyqubitoperators!(ψ, O_, qubitapplymode; N=N, n=n, m=m)

    # LAST STEP: exp(𝒊 HD t[numsteps])), ie. exp(-𝒊 HD T)
    ψ .= UD' * ψ                        # ROTATE INTO DEVICE BASIS
    ψ .*= exp.( (im*T) * ΛD)            # ROTATE PHASES FOR ONE LAST TIME EVOLUTION

    ######################################################################################

    # RE-NORMALIZE THIS STATE
    ψ .= ψ / √abs(ψ'*ψ)

    if iobasis isa QubitBasis;  ψ .= UD  * ψ;   end;    # ROTATE *OUT* OF DEVICE BASIS
end


"""
    _preparequbitdrives(
        pulses::AbstractVector{<:Pulses.PulseTemplate}, # PULSE TEMPLATES FOR EACH QUBIT
        m::Integer,                                     # NUMBER OF LEVELS ON EACH QUBIT
        t::Number,                                      # TIME POINT
        Δt::Number;                                     # TIME TO THE NEXT TIME POINT

        # INFERRED VALUES (relatively fast, but you can pass them in if you'd like)
        n = length(pulses),                             # NUMBER OF QUBITS

        # CALCULATED VALUES (pass these in to significantly speed up optimizations)
        a = Utils.a_matrix(m),                          # SINGLE-QUBIT ANNIHILATION OPERATOR

        # PRE-ALLOCATIONS (for those that want every last drop of efficiency...)
        O_ = [Matrix{ComplexF64}(undef, m,m) for q ∈ 1:n]   # HOLDS QUBIT-WISE DRIVE OPERATORS
    )

Prepare a vector of qubit operations representing the instantaneous action of a pulse.

Say a pulse has amplitude ``Ω_q`` and frequency ``ν_q``,
    and define ``z_q = Ω_q \\exp(𝒊·ν_q·t)``.
We may model the action of the pulse on a resonant system at time ``t``
    with a "Control" Hamiltonian ``H_q = z_q a_q + z̄_q a^\\dagger``,
    and the evolution over short time ``Δt`` as ``\\exp(-𝒊·Δt·H)``.

"""
function _preparequbitdrives(
    pulses::AbstractVector{<:Pulses.PulseTemplate}, # PULSE TEMPLATES FOR EACH QUBIT
    m::Integer,                                     # NUMBER OF LEVELS ON EACH QUBIT
    t::Number,                                      # TIME POINT
    Δt::Number;                                     # TIME TO THE NEXT TIME POINT

    # INFERRED VALUES (relatively fast, but you can pass them in if you'd like)
    n = length(pulses),                             # NUMBER OF QUBITS

    # CALCULATED VALUES (pass these in to significantly speed up optimizations)
    a = nothing,                                    # SINGLE-QUBIT ANNIHILATION OPERATOR

    # PRE-ALLOCATIONS (for those that want every last drop of efficiency...)
    O_ = [Matrix{ComplexF64}(undef, m,m) for q ∈ 1:n]   # HOLDS QUBIT-WISE DRIVE OPERATORS
)
    ######################################################################################
    #                            PRELIMINARY CALCULATIONS

    if a === nothing
        a = Utils.a_matrix(m)                       # SINGLE-QUBIT ANNIHILATION OPERATOR
    end

    ######################################################################################
    #                               PREPARE QUBIT DRIVES
    for q ∈ 1:n
        # EXTRACT TIME-DEPENDENT COEFFICIENTS
        Ω = Pulses.amplitude(pulses[q], t)
        ν = Pulses.frequency(pulses[q], t)
        z = Ω * exp(im*ν*t)

        # CONSTRUCT AND EXPONENTIATE MATRIX
        O_[q] .= exp((-im*Δt) * (z*a + z'*a'))
    end
    ######################################################################################

    return O_
end






"""
    evolve!(
        ψ::AbstractVector{<:Number},
        pulses::AbstractVector{<:Pulses.PulseTemplate},
        device::Devices.Device,
        ::Type{Prediag};
        iobasis::IOBasisMode = DeviceBasis(),
        numsteps::Integer = 2000,
        qubitapplymode::QubitApplyMode = Kronec(),
        suzukiorder::Integer = 2,

        # INFERRED VALUES (relatively fast, but pass them in to minimize allocations)
        N = length(ψ),                      # SIZE OF STATEVECTOR
        n = length(device),                 # NUMBER OF QUBITS
        m = round(Int, N^(1/n)),            # NUMBER OF LEVELS ON EACH QUBIT
        T = length(pulses[1]),              # TOTAL DURATION OF EVOLUTION
        t_= range(0,T,numsteps+1),          # TIME GRID
        Δt= T / numsteps,                   # DURATION OF EACH TIME STEP

        # CALCULATED VALUES (pass these in to significantly speed up optimizations)
        ΛD = nothing,                       # EIGENVALUES OF STATIC HAMILTONIAN
        UD = nothing,                       # CORRESPONDING EIGENVECTORS
        Λ  = nothing,                       # EIGENVALUES OF Q (OR P!) OPERATOR
        UQP= nothing,                       # ROTATION FROM P->Q BASIS
        UPQ= nothing,                       # ROTATION FROM Q->P BASIS
        in_basis = nothing,                 # STARTING BASIS FOR DRIVE OPERATION
        outbasis = nothing,                 #   ENDING BASIS FOR DRIVE OPERATION
        L  = nothing,                       # LIGAND (STATIC PROPAGATION) OPERATION

        # PRE-ALLOCATIONS (for those that want every last drop of efficiency...)
        O_ = [Matrix{ComplexF64}(undef, m,m) for q ∈ 1:n]   # QUBIT-WISE DRIVE OPERATORS
    )

Trotterize the time evolution operator,
    but switch back and forth between the static and drive bases
    so that the time evolution in each is more efficient.
Additionally, decompose the drive Hamiltonian into two time-independent components,
    so that all time evolution can be computed by only
    matrix-vector multiplications and vector exponentiations.

This method is almost identical to the `Rotate` method,
    but evades *all* eigenvalue calculations within a time step.
Thus, it invokes only matrix-vector multiplications (``O(N^2)``)
    and vector exponentiations (``O(N)``) at each time step.
It sounds at first impression like it should be even faster than `Rotate`,
    but it turns out ``N`` is usually bigger than ``n·m^3`` so it's not...
Additionally, it incurs some extra error from the drive Hamiltonian factorization.

The keyword argument `numsteps` specifies the number of time steps;
    time scales linearly, and accuracy scales inversely.

The keyword argument `qubitapplymode` specifies the numerical approach to
    applying the qubit-wise drive operators at each time step.
The choices are:
- `Kronec()` this method combines all qubit-wise operators into a single N×N matrix,
    then applies them with a single matrix-vector multiplication.
- `Tensor()` this method reshapes the statevector into an n-dimensional array,
    and performs a tensor contraction over each qubit-wise operator.
  In principle this one should scale significantly better than `Kronec`,
    but in practice the overhead from tensor manipulation may be steep.
TODO: The horrible overhead on tensor algebra might just vanish by caching...

The keyword argument `suzukiorder` controls the product formula
    used to expand the drive Hamiltonian.
- Suzuki order 1 corresponds to ``\\exp(A·B)≈\\exp(A)·\\exp(B)``
- Suzuki order 2 corresponds to ``\\exp(A·B)≈\\exp(A/2)·\\exp(B)·\\exp(A/2)``
- Suzuki order 0 isn't really a thing but this method uses it to correspond to
    ``\\exp(A·B)≈\\exp(A)·\\exp(B)·\\exp(-[A,B]/2)``.
  In fact, for the decomposition used, this formula is algebraically exact.
  Unfortunately, this exactness vanishes completely when representing
    bosonic algebra with finite matrices, so it's only useful for large `m`,
    which are, alas, computationally inaccessible.
  In other words, don't bother with `suzukiorder=0`...
- Higher-order product formulae are mathematically defined, but not implemented.

"""
function evolve!(
    ψ::AbstractVector{<:Number},
    pulses::AbstractVector{<:Pulses.PulseTemplate},
    device::Devices.Device,
    ::Type{Prediag};
    iobasis::IOBasisMode = DeviceBasis(),
    numsteps::Integer = 2000,
    qubitapplymode::QubitApplyMode = Kronec(),
    suzukiorder::Integer = 2,

    # INFERRED VALUES (relatively fast, but pass them in to minimize allocations)
    N = length(ψ),                      # SIZE OF STATEVECTOR
    n = length(device),                 # NUMBER OF QUBITS
    m = round(Int, N^(1/n)),            # NUMBER OF LEVELS ON EACH QUBIT
    T = length(pulses[1]),              # TOTAL DURATION OF EVOLUTION
    t_= range(0,T,numsteps+1),          # TIME GRID
    Δt= T / numsteps,                   # DURATION OF EACH TIME STEP

    # CALCULATED VALUES (pass these in to significantly speed up optimizations)
    ΛD = nothing,                       # EIGENVALUES OF STATIC HAMILTONIAN
    UD = nothing,                       # CORRESPONDING EIGENVECTORS
    Λ  = nothing,                       # EIGENVALUES OF Q (OR P!) OPERATOR
    UQP= nothing,                       # ROTATION FROM P->Q BASIS
    UPQ= nothing,                       # ROTATION FROM Q->P BASIS
    in_basis = nothing,                 # STARTING BASIS FOR DRIVE OPERATION
    outbasis = nothing,                 #   ENDING BASIS FOR DRIVE OPERATION
    L  = nothing,                       # LIGAND (STATIC PROPAGATION) OPERATION

    # PRE-ALLOCATIONS (for those that want every last drop of efficiency...)
    O_ = [Matrix{ComplexF64}(undef, m,m) for q ∈ 1:n]   # QUBIT-WISE DRIVE OPERATORS
)
    ######################################################################################
    #                            PRELIMINARY CALCULATIONS

    if any((ΛD, UD) .=== nothing)
        HD = Devices.static_hamiltonian(device, m)  # DEVICE HAMILTONIAN
        ΛD, UD = Utils.dressedbasis(HD)             # DEVICE EIGENVALUES AND EIGENVECTORS
    end; if any((Λ, UQP) .=== nothing) || (suzukiorder==2 && UPQ === nothing)
        a = Utils.a_matrix(m)                       # SINGLE-QUBIT ANNIHILATION OPERATOR
        Q =      (a + a')                           # CANONICAL COORDINATE OPERATOR
        P = im * (a - a')                           # CANONICAL   MOMENTUM OPERATOR
        ΛUQ = eigen(Hermitian(Q))
        ΛUP = eigen(Hermitian(P))
        Λ = ΛUQ.values                              # EIGENVALUES OF Q OPERATOR (OR P)
        UQ= ΛUQ.vectors
        UP= ΛUP.vectors
        UQP = UQ' * UP                              # ROTATION MATRIX FROM P -> Q BASIS
        UPQ = UP' * UQ                              # ROTATION MATRIX FROM Q -> P BASIS
    end; if in_basis === nothing
        in_basis = Utils.kron_concat(suzukiorder==2 ? UQ : UP, n)   # DRIVE'S  INPUT BASIS
    end; if outbasis === nothing
        outbasis = Utils.kron_concat(UQ, n)                         # DRIVE'S OUTPUT BASIS
    end; if L === nothing
        L = in_basis' * UD * Diagonal(exp.((-im*Δt)*ΛD)) * UD' * outbasis   # LIGAND OP.
    end

    if iobasis isa DeviceBasis; ψ .= UD * ψ;    end;    # ROTATE *OUT* OF DEVICE BASIS

    ######################################################################################
    #                                 TIME EVOLUTION

    #= NOTE: The very first step is, mathematically, exp(-𝒊 HD t_[1]),
        but since t_[1]=0, this is an identity operation and we can skip it. =#

    # ROTATE INTO `in_basis`
    ψ .= in_basis' * ψ

    # APPLY FIRST PULSE DRIVES  (use Δt/2 for first and last time step)
    _preparequbitdrives_productformula(pulses, m, t_[1], Δt/2; suzukiorder=suzukiorder,
        Λ=Λ, UQP=UQP, UPQ=UPQ, n=n, O_=O_
    )
    _applyqubitoperators!(ψ, O_, qubitapplymode; N=N, n=n, m=m)

    for i ∈ 2:numsteps
        # CONNECT EACH TIME STEP WITH THE DEVICE ACTION
        ψ .= L * ψ

        # APPLY PULSE DRIVES
        _preparequbitdrives_productformula(pulses, m, t_[i], Δt; suzukiorder=suzukiorder,
            Λ=Λ, UQP=UQP, UPQ=UPQ, n=n, O_=O_
        )
        _applyqubitoperators!(ψ, O_, qubitapplymode; N=N, n=n, m=m)
    end

    # APPLY LAST PULSE DRIVES   (use Δt/2 for first and last time step)
    ψ .= L * ψ
    _preparequbitdrives_productformula(pulses, m, t_[end], Δt/2; suzukiorder=suzukiorder,
        Λ=Λ, UQP=UQP, UPQ=UPQ, n=n, O_=O_
    )
    _applyqubitoperators!(ψ, O_, qubitapplymode; N=N, n=n, m=m)

    # ROTATE *OUT* OF `outbasis`
    ψ .= outbasis * ψ

    # LAST STEP: exp(𝒊 HD t[numsteps])), ie. exp(-𝒊 HD T)
    ψ .= UD' * ψ                        # ROTATE INTO DEVICE BASIS
    ψ .*= exp.( (im*T) * ΛD)            # ROTATE PHASES FOR ONE LAST TIME EVOLUTION

    ######################################################################################

    # RE-NORMALIZE THIS STATE
    ψ .= ψ / √abs(ψ'*ψ)

    if iobasis isa QubitBasis;  ψ .= UD  * ψ;   end;    # ROTATE *OUT* OF DEVICE BASIS
end

"""
    _preparequbitdrives_productformula!(
        pulses::AbstractVector{<:Pulses.PulseTemplate}, # PULSE TEMPLATES FOR EACH QUBIT
        m::Integer,                                     # NUMBER OF LEVELS ON EACH QUBIT
        t::Number,                                      # TIME POINT
        Δt::Number;                                     # TIME TO THE NEXT TIME POINT
        suzukiorder = 2,                                # SUZUKI ORDER OF PRODUCT FORMULA

        # MANDATORY (*could* be calculated, but it doesn't seem worth the trouble...)
        Λ   = nothing,                                  # EIGENVALUES OF Q OPERATOR (OR P)
        UQP = nothing,                                  # ROTATION MATRIX FROM P -> Q BASIS
        UPQ = nothing,                                  # ROTATION MATRIX FROM Q -> P BASIS

        # INFERRED VALUES (relatively fast, but you can pass them in if you'd like)
        n = length(pulses),                             # NUMBER OF QUBITS

        # PRE-ALLOCATIONS (for those that want every last drop of efficiency...)
        O_ = [Matrix{ComplexF64}(undef, m,m) for q ∈ 1:n]   # HOLDS QUBIT-WISE DRIVE OPERATORS
    )

Prepare a vector of qubit operations representing the instantaneous action of a pulse.

Say a pulse has amplitude ``Ω_q`` and frequency ``ν_q``,
    and define ``z_q = Ω_q \\exp(𝒊·ν_q·t)``.
We may model the action of the pulse on a resonant system at time ``t``
    with a "Control" Hamiltonian ``H_q = z_q a_q + z̄_q a^\\dagger``,
    and the evolution over short time ``Δt`` as ``\\exp(-𝒊·Δt·H)``.

This variant rewrites the drive Hamiltonian ``H_q → x·Q + y·P``,
    to rewrite the evolution operator ``\\exp(-𝒊·Δt·H)≈\\exp(-𝒊·Δt·x·Q)\\exp(-𝒊·Δt·x·P)``
        or a related product formula, selected with `suzukiorder`.

As an EXTRA feature, `suzukiorder=0` will do a first-order product formula,
    but include the commutator ``\\exp(Δt²·x·y·[Q,P]/2)``, where ``[Q,P]=-2𝒊``.
    In the limit where ``m→∞``, this is exact.
But, uh, we're not in that limit, so...it's just for fun... ^_^

"""
function _preparequbitdrives_productformula(
    pulses::AbstractVector{<:Pulses.PulseTemplate}, # PULSE TEMPLATES FOR EACH QUBIT
    m::Integer,                                     # NUMBER OF LEVELS ON EACH QUBIT
    t::Number,                                      # TIME POINT
    Δt::Number;                                     # TIME TO THE NEXT TIME POINT
    suzukiorder = 2,                                # SUZUKI ORDER OF PRODUCT FORMULA

    # INFERRED VALUES (relatively fast, but pass them in to minimize allocations)
    n = length(pulses),                             # NUMBER OF QUBITS

    # CALCULATED VALUES (pass these in to significantly speed up optimizations)
    Λ   = nothing,                                  # EIGENVALUES OF Q OPERATOR (OR P)
    UQP = nothing,                                  # ROTATION MATRIX FROM P -> Q BASIS
    UPQ = nothing,                                  # ROTATION MATRIX FROM Q -> P BASIS

    # PRE-ALLOCATIONS (for those that want every last drop of efficiency...)
    O_ = [Matrix{ComplexF64}(undef, m,m) for q ∈ 1:n]   # HOLDS QUBIT-WISE DRIVE OPERATORS
)
    ######################################################################################
    #                            PRELIMINARY CALCULATIONS

    if any((Λ, UQP) .=== nothing) || (suzukiorder==2 && UPQ === nothing)
        a = Utils.a_matrix(m)                       # SINGLE-QUBIT ANNIHILATION OPERATOR
        Q =      (a + a')                           # CANONICAL COORDINATE OPERATOR
        P = im * (a - a')                           # CANONICAL   MOMENTUM OPERATOR
        ΛUQ = eigen(Hermitian(Q))
        ΛUP = eigen(Hermitian(P))
        Λ = ΛUQ.values                              # EIGENVALUES OF Q OPERATOR (OR P)
        UQ= ΛUQ.vectors
        UP= ΛUP.vectors
        UQP = UQ' * UP                              # ROTATION MATRIX FROM P -> Q BASIS
        UPQ = UP' * UQ                              # ROTATION MATRIX FROM Q -> P BASIS
    end;

    ######################################################################################
    #                               PREPARE QUBIT DRIVES

    for q ∈ 1:n
        # EXTRACT TIME-DEPENDENT COEFFICIENTS
        Ω = Pulses.amplitude(pulses[q], t)
        ν = Pulses.frequency(pulses[q], t)
        z = Ω * exp(im*ν*t)
        x, y = real(z), imag(z)

        # EVOLVE QUBIT IN TIME, AND EXTEND FULL-QUBIT OPERATOR
        if     suzukiorder == 0
            expQ = Diagonal(exp.((-im*Δt*x) * Λ))
            expP = Diagonal(exp.((-im*Δt*y) * Λ))

            O_[q] .= expQ * UQP * expP * exp(-im*x*y*Δt^2)
                # Alas, this is only going to work for large m.
        elseif suzukiorder == 1
            expQ = Diagonal(exp.((-im*Δt*x) * Λ))
            expP = Diagonal(exp.((-im*Δt*y) * Λ))

            O_[q] .= expQ * UQP * expP
        elseif suzukiorder == 2
            expQ = Diagonal(exp.((-im*Δt*x/2) * Λ))
            expP = Diagonal(exp.((-im*Δt*y  ) * Λ))

            O_[q] .= expQ * UQP * expP * UPQ * expQ
        else
            error("Only `suzukiorder`s 0, 1, and 2 are supported.")
        end
    end

    # TODO: Pre-allocate expQ, expP

    ######################################################################################
    return O_
end




"""
    _applyqubitoperators!(
        ψ::AbstractVector{<:Number},
        O_::AbstractVector{<:AbstractMatrix{<:Number}},
        mode::QubitApplyMode;

        # INFERRED VALUES (relatively fast, but you can pass them in if you'd like)
        #= None of these are needed for mode::Kronec, you may pass them in as `nothing` =#
        N = length(ψ),                          # SIZE OF STATEVECTOR
        n = length(O_),                         # NUMBER OF QUBITS
        m = round(Int, N^(1/n)),                # NUMBER OF LEVELS ON EACH QUBIT
    )

Apply a sequence of qubit operators to a statevector ψ.

In principle, the most efficient means to do this is to reshape ψ as an n-body tensor,
    and opply each operator with a tensor contraction.
To do this, use `mode=Tensor()`

But there's a significant amount of overhead,
    (most likely due to copies made during permutation of dimensions),
    so in many (if not most) cases, it is preferable to use `mode=Kronec()`.
This will simply combine all qubit operators into a full N×N matrix and left-multiply `ψ`.

TODO: The horrible overhead on tensor algebra might just vanish by caching...

"""
function _applyqubitoperators!(
    ψ::AbstractVector{<:Number},
    O_::AbstractVector{<:AbstractMatrix{<:Number}},
    mode::QubitApplyMode;

    # INFERRED VALUES (relatively fast, but pass them in to minimize allocations)
    N = length(ψ),                          # SIZE OF STATEVECTOR
    n = length(O_),                         # NUMBER OF QUBITS
    m = round(Int, N^(1/n)),                # NUMBER OF LEVELS ON EACH QUBIT
)
    ######################################################################################
    if     mode isa Kronec
        O = Utils.kron_concat(O_)               # FULL HILBERT-SPACE OPERATOR
        ψ .= O * ψ                              # APPLIED TO STATEVECTOR
        # TODO: Pre-allocate O..?

    ######################################################################################
    elseif mode isa Tensor

        ψ_ = reshape(ψ, (m for _ in 1:n)...)    # *NOT* A COPY; MUTATIONS APPLY TO BOTH
        ψ_ .= ncon(
            [O_..., ψ_],                            # LIST OF TENSORS
            [([-q, q] for q in 1:n)..., n:-1:1],    # LIST OF INDICES ON EACH TENSOR
            # zeros(Bool, n+1), :cache,               # ENABLE CACHING
            output=-n:-1,                           # FINAL PERMUTATION
        )
        # ψ HAS ALREADY BEEN UPDATED, IN MUTATIONS OF ψ_

        #= TODO: Proper benchmarking to understand if caching is useful.
            It obviously helps over repeated trials in @btime; that's not fair.
            It *probably* helps over repeated contractions over many time steps,
                so it's surely worth having.
            But I'd like to understand how it works better before I permit it.
                Does it speed allocations *within* a single contraction?
        =#

    ######################################################################################
    else
        error("Invalid `QubitApplyMode` object. (How did you manage that???)")
    end
end

end # END MODULE
