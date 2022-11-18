#= Code to evolve a quantum-controlled system in time. =#

module Evolutions

import LinearAlgebra: eigen, Hermitian, Diagonal, norm, lmul!, rmul!, mul!
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

#= TODO: There are a couple points where we were stubbornly resisting
        the "flexible methods, rigid structs" paradigm.
    Double check that looser methods don't slow you down, then go all in.
=#

#= TODO: Structs should allow for alternate float precision. =#

#= TODO: Match notation in notebook. HIC→V, V→L, expHIC→E, etc.

Also temp variables look awkward.
Do something like `_` prefix to indicate it's an allocation,
    then `1D` or `2D` for the vector vs. matrix distinction,
    then `N` or `n` to indicate size,
    then finally a `_` suffix to indicate array.
`tmpK_` can change to something more innocuous like `_p_`
=#

#= TODO: I/O and QubitApply can be :symbols instead of types, I think. =#

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
        ψ::Vector{ComplexF64},
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
    ψ::Vector{ComplexF64},
    pulses::AbstractVector{<:Pulses.PulseTemplate},
    device::Devices.Device;
    iobasis::IOBasisMode = DeviceBasis()
) = evolve!(ψ, pulses, device, Rotate; iobasis=iobasis)




"""
    evolve!(
        ψ::Vector{ComplexF64},
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

        # CALCULATED VALUES (pass these in to significantly speed up optimizations)
        tmpM = Matrix{ComplexF64}(undef, N,N),          # FOR CONTROL HAMILTONIAN
        tmpV = Vector{ComplexF64}(undef, N),            # FOR DEVICE ACTION
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
    ψ::Vector{ComplexF64},
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

    # PRE-ALLOCATIONS (for those that want every last drop of efficiency...)
    tmpM = Matrix{ComplexF64}(undef, N,N),          # FOR CONTROL HAMILTONIAN
    tmpV = Vector{ComplexF64}(undef, N),            # FOR DEVICE ACTION
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

    # ROTATE INTO DEVICE BASIS
    if iobasis isa QubitBasis;  Utils.transform!(ψ, UD', tmpV); end

    ######################################################################################
    #                       DEFINE AND SOLVE DIFFERENTIAL EQUATIONS

    # SOLVE THE SYSTEM
    #= NOTE:
        This method autoselects an algorithm based on I have no idea what,
            meaning I have no idea what the time complexity or accuracy are likely to be,
            or how I should expect them to scale with increasing system size.
        But, it works *pretty* well for the single-qubit case,
            so I'm willing to treat it as a sort of black-box standard.
    =#
    p = (pulses, N, n, ΛD, a_, tmpM, tmpV)
    schrodinger = ODEProblem(_interaction!, ψ, (0.0, T), p)
    solution = solve(schrodinger, save_everystep=false)     # TIME-CONSUMING STEP

    # WRITE FINAL SOLUTION TO THE GIVEN STATEVECTOR
    ψ .= solution.u[end]

    # RE-NORMALIZE THIS STATE
    ψ ./= norm(ψ)

    ######################################################################################

    # ROTATE *OUT* OF DEVICE BASIS
    if iobasis isa QubitBasis;  Utils.transform!(ψ, UD, tmpV)   end
end

""" Auxiliary function for Direct `evolve!`.
Mutates derivative vector du with Schrodinger's interaction picture equation.
p is a tuple of "parameters", actually just used to pass constants/preallocations.
"""
function _interaction!(du, u, p, t)
    # UNPACK PARAMETERS
    pulses  = p[1]
    N       = p[2]
    n       = p[3]
    ΛD      = p[4]
    a_      = p[5]
    tmpM    = p[6]
    tmpV    = p[7]

    # CONSTRUCT CONTROL HAMILTONIAN (IN DEVICE BASIS)
    tmpM .= zeros(N,N)
    for q ∈ 1:n
        Ω = Pulses.amplitude(pulses[q], t)
        ν = Pulses.frequency(pulses[q], t)
        z = Ω * exp(im*ν*t)

        tmpM .+= z .* a_[q]     # ADD IN za terms
    end
    tmpM .+= tmpM'              # ADD IN z*a† terms

    # CONSTRUCT INTERACTION PICTURE HAMILTONIAN
    tmpV .= exp.((im*t) .* ΛD)                  # DEVICE ACTION
    expD = Diagonal(tmpV)                           # CREATE A DIAGONAL-MATRIX VIEW
    lmul!(expD, tmpM); rmul!(tmpM, expD')       # CONJUGATE WITH DEVICE ACTION

    # SCHRODINGER'S EQUATION
    lmul!(-im, tmpM)                            # ADD THE i FROM SCHRODINGER'S EQN
    mul!(du, tmpM, u)                           # SET du = -i H_I u
end



"""
    evolve!(
        ψ::Vector{ComplexF64},
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

        # PRE-ALLOCATIONS (for those that want every last drop of efficiency...)
        tmpM = Matrix{ComplexF64}(undef, N,N),          # FOR CONTROL HAMILTONIAN
        tmpV = Vector{ComplexF64}(undef, N),            # FOR DEVICE ACTION
    )

Trotterize the time-evolution operator,
    directly exponentiating ``\\exp(-𝒊·Δt·H)`` at each time step.

The keyword argument `numsteps` specifies the number of time steps;
    time scales linearly, and accuracy scales inversely.

This method is the only one that invokes an ``O(N^3)`` operation at every time step.
Don't use it, except to illustrate how much better other methods are. ^_^

"""
function evolve!(
    ψ::Vector{ComplexF64},
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

    # PRE-ALLOCATIONS (for those that want every last drop of efficiency...)
    tmpM = Matrix{ComplexF64}(undef, N,N),          # FOR CONTROL HAMILTONIAN
    tmpV = Vector{ComplexF64}(undef, N),            # FOR DEVICE ACTION
)
    ######################################################################################
    #                            PRELIMINARY CALCULATIONS

    if any((ΛD, UD) .=== nothing)
        HD = Devices.static_hamiltonian(device, m)  # DEVICE HAMILTONIAN
        ΛD, UD = Utils.dressedbasis(HD)             # DEVICE EIGENVALUES AND EIGENVECTORS
    end; if a_ === nothing
        a_ = Utils.algebra(n, m, basis=UD)          # LIST OF ROTATED ANNIHILATION OPERATORS
    end

    # ROTATE INTO DEVICE BASIS
    if iobasis isa QubitBasis;  Utils.transform!(ψ, UD', tmpV); end

    ######################################################################################
    #                                 TIME EVOLUTION

    # FIRST TIME STEP   (use Δt/2 for first and last time step)
    ψ .= _step(ψ, t_[1], Δt/2, Direct, pulses, N, n, ΛD, a_, tmpM, tmpV)

    for i ∈ 2:numsteps
        ψ .= _step(ψ, t_[i], Δt, Direct, pulses, N, n, ΛD, a_, tmpM, tmpV)
    end

    # LAST TIME STEP    (use Δt/2 for first and last time step)
    ψ .= _step(ψ, t_[end], Δt/2, Direct, pulses, N, n, ΛD, a_, tmpM, tmpV)

    ######################################################################################

    # RE-NORMALIZE THIS STATE
    ψ ./= norm(ψ)

    # ROTATE *OUT* OF DEVICE BASIS
    if iobasis isa QubitBasis;  Utils.transform!(ψ, UD, tmpV)   end
end

""" Auxiliary function for Direct `evolve!`. """
function _step(ψ, t, Δt, ::Type{Direct}, pulses, N, n, ΛD, a_, tmpM, tmpV)
    ######################################################################################
    #                                 SINGLE TIME STEP

    # CONSTRUCT CONTROL HAMILTONIAN (IN DEVICE BASIS)
    tmpM .= zeros(N,N)
    for q ∈ 1:n
        Ω = Pulses.amplitude(pulses[q], t)
        ν = Pulses.frequency(pulses[q], t)
        z = Ω * exp(im*ν*t)

        tmpM .+= z .* a_[q]     # ADD IN za terms
    end
    tmpM .+= tmpM'              # ADD IN z*a† terms

    # CONSTRUCT INTERACTION PICTURE HAMILTONIAN
    tmpV .= exp.((im*t) .* ΛD)                  # DEVICE ACTION
    expD = Diagonal(tmpV)                           # CREATE A DIAGONAL-MATRIX VIEW
    lmul!(expD, tmpM); rmul!(tmpM, expD')       # CONJUGATE WITH DEVICE ACTION

    # PREPARE TIME-EVOLUTION OPERATOR
    tmpM .= exp((-im*Δt).*tmpM)
    #= NOTE: *THIS* step is the bottleneck,
                not only because it is algorithmically the most complex,
                but because it is the *ONLY* step inside the time loop
                that invokes any allocations at all.
                Interestingly, I think I can cut down those allocations by a factor of two
                by manually diagonalizing `tmpM`
                then doing in-place matrix multiplications.

                In fact, `exp` seems to be doing a different algorithm entirely,
                since...its argument is not _Hermitian_.
            ...it's *anti*-Hermitian! That is dumb. Maybe worth filing an issue over.
                Anyways, the algorithm being used is called "squaring and scaling".
                I dunno how it works but it seems to be a bit slower,
                at least for the limited N=100 test I did in the REPL.
                The point is we can probably cut time in two by manually exponentiating.
    =#

    # APPLY TIME-EVOLUTION OPERATOR
    mul!(tmpV, tmpM, ψ)
    return tmpV
end




"""
    evolve!(
        ψ::Vector{ComplexF64},
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

        # PRE-ALLOCATIONS (for those that want every last drop of efficiency...)
        tmpM = Matrix{ComplexF64}(undef, N,N),          # FOR CONTROL HAMILTONIAN
        tmpV = Vector{ComplexF64}(undef, N),            # FOR DEVICE ACTION
    )

Trotterize the time-evolution operator,
    calculating the matrix exponential action ``\\exp(-𝒊·Δt·H) |ψ⟩`` at each time step.

The keyword argument `numsteps` specifies the number of time steps;
    time scales linearly, and accuracy scales inversely.

This method is the only one that invokes an ``O(N^3)`` operation at every time step.
Don't use it, except to illustrate how much better other methods are. ^_^

"""
function evolve!(
    ψ::Vector{ComplexF64},
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

    # PRE-ALLOCATIONS (for those that want every last drop of efficiency...)
    tmpM = Matrix{ComplexF64}(undef, N,N),          # FOR CONTROL HAMILTONIAN
    tmpV = Vector{ComplexF64}(undef, N),            # FOR DEVICE ACTION
)
    ######################################################################################
    #                            PRELIMINARY CALCULATIONS

    if any((ΛD, UD) .=== nothing)
        HD = Devices.static_hamiltonian(device, m)  # DEVICE HAMILTONIAN
        ΛD, UD = Utils.dressedbasis(HD)             # DEVICE EIGENVALUES AND EIGENVECTORS
    end; if a_ === nothing
        a_ = Utils.algebra(n, m, basis=UD)          # LIST OF ROTATED ANNIHILATION OPERATORS
    end

    # ROTATE INTO DEVICE BASIS
    if iobasis isa QubitBasis;  Utils.transform!(ψ, UD', tmpV); end

    ######################################################################################
    #                                 TIME EVOLUTION

    # FIRST TIME STEP   (use Δt/2 for first and last time step)
    ψ .= _step(ψ, t_[1], Δt/2, Lanczos, pulses, N, n, ΛD, a_, tmpM, tmpV)

    for i ∈ 2:numsteps
        ψ .= _step(ψ, t_[i], Δt, Lanczos, pulses, N, n, ΛD, a_, tmpM, tmpV)
    end

    # LAST TIME STEP    (use Δt/2 for first and last time step)
    ψ .= _step(ψ, t_[end], Δt/2, Lanczos, pulses, N, n, ΛD, a_, tmpM, tmpV)

    ######################################################################################

    # RE-NORMALIZE THIS STATE
    ψ ./= norm(ψ)

    # ROTATE *OUT* OF DEVICE BASIS
    if iobasis isa QubitBasis;  Utils.transform!(ψ, UD, tmpV);  end
end


""" Auxiliary function for Lanczos `evolve!`. """
function _step(ψ, t, Δt, ::Type{Lanczos}, pulses, N, n, ΛD, a_, tmpM, tmpV)
    ######################################################################################
    #                                 SINGLE TIME STEP

    # CONSTRUCT CONTROL HAMILTONIAN (IN DEVICE BASIS)
    tmpM .= zeros(N,N)
    for q ∈ 1:n
        Ω = Pulses.amplitude(pulses[q], t)
        ν = Pulses.frequency(pulses[q], t)
        z = Ω * exp(im*ν*t)

        tmpM .+= z .* a_[q]     # ADD IN za terms
    end
    tmpM .+= tmpM'              # ADD IN z*a† terms

    # CONSTRUCT INTERACTION PICTURE HAMILTONIAN
    tmpV .= exp.((im*t) .* ΛD)                  # DEVICE ACTION
    expD = Diagonal(tmpV)                           # CREATE A DIAGONAL-MATRIX VIEW
    lmul!(expD, tmpM); rmul!(tmpM, expD')       # CONJUGATE WITH DEVICE ACTION

    # APPLY TIME-EVOLUTION OPERATOR
    tmpV = exponentiate(tmpM, -im*Δt, ψ)[1]
    #= NOTE: *THIS* step is the bottleneck,
                not only because it is algorithmically the most complex,
                but because it is the *ONLY* step inside the time loop
                that invokes any allocations at all.
    =#

    return tmpV
end




"""
    evolve!(
        ψ::Vector{ComplexF64},
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
        tmpV = Vector{ComplexF64}(undef, N),                # FOR MATRIX-VECTOR MULTIPLICATION
        tmpM_ = [Matrix{ComplexF64}(undef, m,m) for q ∈ 1:n],   # QUBIT-WISE DRIVE OPERATORS
        tmpK_ = nothing,                                    # FOR APPLYING OPERATORS
                                                        # (default depends on `qubitapplymode`)
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

"""
function evolve!(
    ψ::Vector{ComplexF64},
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
    tmpV = Vector{ComplexF64}(undef, N),                # FOR MATRIX-VECTOR MULTIPLICATION
    tmpM_ = [Matrix{ComplexF64}(undef, m,m) for q ∈ 1:n],   # QUBIT-WISE DRIVE OPERATORS
    tmpK_ = nothing,                                    # FOR APPLYING OPERATORS
                                                    # (default depends on `qubitapplymode`)
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
    end; if tmpK_ === nothing
        tmpK_ = (
            qubitapplymode isa Kronec ?
                [Matrix{ComplexF64}(undef, m^q, m^q) for q ∈ 1:n] :
            qubitapplymode isa Tensor ? [
                Dims(m for _ in 1:n),                   # RESHAPING DIMENSIONS
                [[[-q, q] for q in 1:n]..., n:-1:1],    # TENSOR INDICES
                -n:-1,                                  # OUTPUT PERMUTATION
                zeros(Bool, n+1),                       # ADJOINT FLAG
            ] : error("Invalid `QubitApplyMode` object. (How did you manage that???)")
        )
    end

    # ROTATE *OUT* OF DEVICE BASIS
    if iobasis isa DeviceBasis; Utils.transform!(ψ, UD, tmpV);  end

    ######################################################################################
    #                                 TIME EVOLUTION

    #= NOTE: The very first step is, mathematically, exp(-𝒊 HD t_[1]),
        but since t_[1]=0, this is an identity operation and we can skip it. =#

    # APPLY FIRST QUBIT DRIVES  (use Δt/2 for first and last time step)
    ψ .= _step(ψ, t_[1], Δt/2, Rotate, pulses, qubitapplymode, n, a, tmpV, tmpM_, tmpK_)

    for i ∈ 2:numsteps
        Utils.transform!(ψ, V, tmpV)        # CONNECT QUBIT DRIVES WITH THE DEVICE ACTION
        ψ .= _step(ψ, t_[i], Δt, Rotate, pulses, qubitapplymode, n, a, tmpV, tmpM_, tmpK_)
    end
    Utils.transform!(ψ, V, tmpV)        # CONNECT QUBIT DRIVES WITH THE DEVICE ACTION

    # APPLY LAST PULSE DRIVES   (use Δt/2 for first and last time step)
    ψ .= _step(ψ, t_[end], Δt/2, Rotate, pulses, qubitapplymode, n, a, tmpV, tmpM_, tmpK_)

    # LAST STEP: exp(𝒊 HD t[numsteps])), ie. exp(-𝒊 HD T)
    Utils.transform!(ψ, UD', tmpV)      # ROTATE INTO DEVICE BASIS
    ψ .*= exp.( (im*T) * ΛD)            # ROTATE PHASES FOR ONE LAST TIME EVOLUTION

    ######################################################################################

    # RE-NORMALIZE THIS STATE
    ψ ./= norm(ψ)

    # ROTATE *OUT* OF DEVICE BASIS
    if iobasis isa QubitBasis;  Utils.transform!(ψ, UD, tmpV);  end
end

""" Auxiliary function for Prediag `evolve!`. """
function _step(ψ, t, Δt, ::Type{Rotate}, pulses, qubitapplymode, n, a, tmpV, tmpM_, tmpK_)
    ######################################################################################
    #                                 SINGLE TIME STEP

    # PREPARE QUBIT DRIVES
    for q ∈ 1:n
        # EXTRACT TIME-DEPENDENT COEFFICIENTS
        Ω = Pulses.amplitude(pulses[q], t)
        ν = Pulses.frequency(pulses[q], t)
        z = Ω * exp(im*ν*t)

        # CONSTRUCT AND EXPONENTIATE MATRIX
        tmpM_[q] .= z .* a  # ADD za TERM
                            # THE z' a' TERM IS ACCOUNTED FOR BY THE `Hermitian` VIEW

        tmpM_[q] .= exp((-im*Δt) .* Hermitian(tmpM_[q]))
            # THIS LAST STEP SHOULD BE THE ONLY ONE REQUIRING ANY ALLOCATIONS
    end

    # APPLY QUBIT DRIVES
    if qubitapplymode isa Kronec
        # KRONECKER MODE: CONSTRUCT FULL-BODY OPERATOR
        O = Utils.kron_concat(tmpM_, tmpK_)
        return mul!(tmpV, O, ψ)
    elseif qubitapplymode isa Tensor
        # TENSOR MODE: RESHAPE AND CONTRACT
        #= TODO: Write manual tensor contraction so you can control pre-allocations.
            I still don't really know what cache is doing,
                but it's not doing everything it could.
            Don't forget to change Prediag call also.
        =#
        ψ_ = reshape(ψ, tmpK_[1])   # *NOT* A COPY; MUTATIONS APPLY TO BOTH
        ψ_ .= ncon(
            [tmpM_..., ψ_],                         # LIST OF TENSORS
            tmpK_[2],    # LIST OF INDICES ON EACH TENSOR
            tmpK_[4], :cache,                       # ENABLE CACHING
            output=tmpK_[3],                        # FINAL PERMUTATION
        )
        # ψ HAS ALREADY BEEN UPDATED, IN MUTATIONS OF ψ_
        return ψ
    else
        error("Invalid `QubitApplyMode` object. (How did you manage that???)")
    end
end






"""
    evolve!(
        ψ::Vector{ComplexF64},
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
        tmpV = Vector{ComplexF64}(undef, N),                # FOR MATRIX-VECTOR MULTIPLICATION
        tmpD = Vector{ComplexF64}(undef, m),                # FOR QUBIT DRIVE DIAGONAL FACTORS
        tmpM_ = [Matrix{ComplexF64}(undef, m,m) for q ∈ 1:n],   # QUBIT-WISE DRIVE OPERATORS
        tmpM = nothing,                                     # EXTRA, FOR `suzukiorder=2` ONLY
        tmpK_ = nothing,                                    # FOR APPLYING OPERATORS
                                                        # (default depends on `qubitapplymode`)
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
    ψ::Vector{ComplexF64},
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
    tmpV = Vector{ComplexF64}(undef, N),                # FOR MATRIX-VECTOR MULTIPLICATION
    tmpD = Vector{ComplexF64}(undef, m),                # FOR QUBIT DRIVE DIAGONAL FACTORS
    tmpM_ = [Matrix{ComplexF64}(undef, m,m) for q ∈ 1:n],   # QUBIT-WISE DRIVE OPERATORS
    tmpM = nothing,                                     # EXTRA, FOR `suzukiorder=2` ONLY
    tmpK_ = nothing,                                    # FOR APPLYING OPERATORS
                                                    # (default depends on `qubitapplymode`)
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
    end; if tmpK_ === nothing
        tmpK_ = (
            qubitapplymode isa Kronec ?
                [Matrix{ComplexF64}(undef, m^q, m^q) for q ∈ 1:n] :
            qubitapplymode isa Tensor ? [
                Dims(m for _ in 1:n),                   # RESHAPING DIMENSIONS
                [[[-q, q] for q in 1:n]..., n:-1:1],    # TENSOR INDICES
                -n:-1,                                  # OUTPUT PERMUTATION
                zeros(Bool, n+1),                       # ADJOINT FLAG
            ] : error("Invalid `QubitApplyMode` object. (How did you manage that???)")
        )
    end; if tmpM === nothing && suzukiorder==2
        tmpM = Matrix{ComplexF64}(undef, m,m)
    end;

    # ROTATE *OUT* OF DEVICE BASIS
    if iobasis isa DeviceBasis; Utils.transform!(ψ, UD, tmpV);  end

    ######################################################################################
    #                                 TIME EVOLUTION

    #= NOTE: The very first step is, mathematically, exp(-𝒊 HD t_[1]),
        but since t_[1]=0, this is an identity operation and we can skip it. =#

    # APPLY FIRST QUBIT DRIVES  (use Δt/2 for first and last time step)
    Utils.transform!(ψ, in_basis', tmpV)
    ψ .= _step(ψ, t_[1], Δt/2, Prediag,
        pulses, suzukiorder, qubitapplymode,    # PARAMETERS
        n, Λ, UQP, UPQ,                         # INFERRED
        tmpV, tmpD, tmpM_, tmpM, tmpK_,         # PRE-ALLOCATION
    )

    for i ∈ 2:numsteps
        Utils.transform!(ψ, L, tmpV)        # CONNECT QUBIT DRIVES WITH THE DEVICE ACTION
        ψ .= _step(ψ, t_[i], Δt, Prediag,       # APPLY QUBIT DRIVES
            pulses, suzukiorder, qubitapplymode,    # PARAMETERS
            n, Λ, UQP, UPQ,                         # INFERRED
            tmpV, tmpD, tmpM_, tmpM, tmpK_,         # PRE-ALLOCATION
        )
    end
    Utils.transform!(ψ, L, tmpV)            # CONNECT QUBIT DRIVES WITH THE DEVICE ACTION

    # APPLY LAST PULSE DRIVES   (use Δt/2 for first and last time step)
    ψ .= _step(ψ, t_[end], Δt/2, Prediag,
        pulses, suzukiorder, qubitapplymode,    # PARAMETERS
        n, Λ, UQP, UPQ,                         # INFERRED
        tmpV, tmpD, tmpM_, tmpM, tmpK_,         # PRE-ALLOCATION
    )
    Utils.transform!(ψ, outbasis, tmpV)

    # LAST STEP: exp(𝒊 HD t[numsteps])), ie. exp(-𝒊 HD T)
    Utils.transform!(ψ, UD', tmpV)      # ROTATE INTO DEVICE BASIS
    ψ .*= exp.( (im*T) .* ΛD)           # ROTATE PHASES FOR ONE LAST TIME EVOLUTION

    ######################################################################################

    # RE-NORMALIZE THIS STATE
    ψ ./= norm(ψ)

    # ROTATE *OUT* OF DEVICE BASIS
    if iobasis isa QubitBasis;  Utils.transform!(ψ, UD, tmpV);  end
end

""" Auxiliary function for Prediag `evolve!`. """
function _step(ψ, t, Δt, ::Type{Prediag},
    pulses, suzukiorder, qubitapplymode,
    n, Λ, UQP, UPQ,
    tmpV, tmpD, tmpM_, tmpM, tmpK_,
)
    ######################################################################################
    #                                 SINGLE TIME STEP

    # PREPARE QUBIT DRIVES
    for q ∈ 1:n
        # EXTRACT TIME-DEPENDENT COEFFICIENTS
        Ω = Pulses.amplitude(pulses[q], t)
        ν = Pulses.frequency(pulses[q], t)
        z = Ω * exp(im*ν*t)
        x, y = real(z), imag(z)

        # COMBINE Q DRIVE, P DRIVE, AND ROTATIONS BETWEEN THEM
        if     suzukiorder <  2
            # CORE ROTATION: Q <- P
            tmpM_[q] .= UQP                 # CORE ROTATION: P -> Q

            # RIGHT-MULTIPLY BY P DRIVE
            tmpD .= exp.((-im*Δt*y) .* Λ)   # P DRIVE CALCULATION
            expP = Diagonal(tmpD)               # DIAGONAL VIEW
            rmul!(tmpM_[q], expP)           # MERGE INTO QUBIT OPERATOR

            #  LEFT-MULTIPLY BY Q DRIVE
            tmpD .= exp.((-im*Δt*x) .* Λ)   # Q DRIVE CALCULATION
            expQ = Diagonal(tmpD)               # DIAGONAL VIEW
            lmul!(expQ, tmpM_[q])           # MERGE INTO QUBIT OPERATOR

            # SPECIAL `suzukiorder=0` MODE: INCLUDE COMMUTATOR
            suzukiorder == 0 && rmul!(tmpM_[q], exp(-im*x*y*Δt^2))
                # Alas, this is only going to work for large m.

        elseif suzukiorder == 2
            # CORE ROTATION: Q <- P
            tmpM .= UQP

            # RIGHT-MULTIPLY BY P DRIVE
            tmpD .= exp.((-im*Δt*y) .* Λ)   # P DRIVE CALCULATION
            expP = Diagonal(tmpD)               # DIAGONAL VIEW
            rmul!(tmpM, expP)               # MERGE INTO QUBIT OPERATOR

            # ADD ROTATION: P <- Q
            mul!(tmpM_[q], tmpM, UPQ)

            # LEFT- AND RIGHT-MULTIPLY BY Q DRIVE (HALF EACH)
            tmpD .= exp.((-im*Δt*x/2) .* Λ) # Q DRIVE CALCULATION
            expQ = Diagonal(tmpD)               # DIAGONAL VIEW
            lmul!(expQ, tmpM_[q]); rmul!(tmpM_[q], expQ)
        else
            error("Only `suzukiorder`s 0, 1, and 2 are supported.")
        end
    end

    # APPLY QUBIT DRIVES
    if qubitapplymode isa Kronec
        # KRONECKER MODE: CONSTRUCT FULL-BODY OPERATOR
        O = Utils.kron_concat(tmpM_, tmpK_)
        return mul!(tmpV, O, ψ)
    elseif qubitapplymode isa Tensor
        # TENSOR MODE: RESHAPE AND CONTRACT
        ψ_ = reshape(ψ, tmpK_[1])   # *NOT* A COPY; MUTATIONS APPLY TO BOTH
        ψ_ .= ncon(
            [tmpM_..., ψ_],                         # LIST OF TENSORS
            tmpK_[2],    # LIST OF INDICES ON EACH TENSOR
            tmpK_[4], :cache,                       # ENABLE CACHING
            output=tmpK_[3],                        # FINAL PERMUTATION
        )
        # ψ HAS ALREADY BEEN UPDATED, IN MUTATIONS OF ψ_
        return ψ
    else
        error("Invalid `QubitApplyMode` object. (How did you manage that???)")
    end
end

end # END MODULE
