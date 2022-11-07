#= Code to evolve a quantum-controlled system in time.

TODO: Just occurred to me, do we need to tell Julia our matrices are Hermitian for eigen()?



For now, treat all input and output vectors as though already in the device basis.
    I'm not convinced this makes sense; I don't think you should need
        or want to conjugate the molecular Hamiltonian.
    The `evolve` function needs to know the U(device<->qubit) matrix *anyways*,
        so it can always *do* the evolution in the device basis
        (preferable perhaps, since in the absence of any pulse, only phases rotate)
        (although it doesn't make any difference at all, in the pre-diagonalized setup).
    But I don't think there's anything *outside* this method to prefer that basis.

    That said, preliminary steps like diagonalizing the device Hamiltonian
        and conjugating the a matrices *should* generally happen *outside* `evolve`,
        since `evolve` will be called many many times for the same device,
            in a typical ctrlVQE optimization.



Evolution Modes:
1. Direct: each time step is independent, calculate U(t->t+Δt) directly
2. Lanczos: avoid exponentiation by calculating action of matrix exponential on vector
3. Rotate: toggle back and forth between device and drive bases
        The device action is the same at each time step,
            and the drive action can be treated in the qubit basis.
4. Prediag: same as rotate, but now apply a product formula on drive Hamiltonian
        This means you no longer need to perform any matrix exponentiation,
            since the eigenvectors of each factor remain static.
        But it does introduce a little extra error, that should vanish as Δt -> 0.

    In both Rotate methods, provide a keyword parameter `tensor` to determine whether
        qubit-wise actions are performed with matrix algebra (low overhead, bad scaling)
        or tensor algebra (optimal scaling but horrible overhead).

        NOTE: the horrible overhead on tensor algebra might just vanish by caching...

=#

# include("./utils.jl")
# include("./pulse.jl")
# include("./device.jl")

module Evolutions
using DifferentialEquations
using LinearAlgebra
using TensorOperations
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

# ENUMERATIONS TO CONTROL QUBIT-WISE OPERATIONS
abstract type QubitApplyMode end
struct Kronec <: QubitApplyMode end
struct Tensor <: QubitApplyMode end

"""
    evolve!(
        ψ::AbstractVector{<:Number},
        pulses::AbstractVector{<:Pulses.PulseTemplate},
        device::Devices.Device;
        numsteps::Integer = 2000
    )

Evolve the state `ψ` in time.

The amount of time evolved is determined by the duration of the pulses,
    which are assumed to have equal duration.

# Arguments
- `ψI` initial statevector of `n>0` qubits each with `nstates` levels
- `pulses` vector of `n` pulse templates
- `device` the `n`-qubit device giving qubit frequencies and couplings
- `numsteps` the number of discrete time units to simulate (ie. Trotter steps)
             must be positive integer

Optionally, specify an EvolutionMode as a final positional argument.
Each mode has a different "setup" phase,
    which can be skipped by passing in keyword arguments.
Some of these setups involve significantly more expensive calculations than
        must occur in the evolution itself,
    so it is highly recommended that you do so.
Study the method headers for each mode in the code for the requisite details.

"""
evolve!(
    ψ::AbstractVector{<:Number},
    pulses::AbstractVector{<:Pulses.PulseTemplate},
    device::Devices.Device;
    numsteps::Integer = 2000
) = evolve!(ψ, pulses, device, Direct; numsteps=numsteps)



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


# TODO: include phase factor from commutator in Prediag mode

# TODO: Treat first and last time-steps differently, using trapezoidal rule.

# TODO: Use the dressed basis scheme from ctrlq.
# TODO: move keyword calculations to `if nothing` lines.




function evolve!(
    ψ::AbstractVector{<:Number},
    pulses::AbstractVector{<:Pulses.PulseTemplate},
    device::Devices.Device,
    ::Type{ODE};

    # INFERRED VALUES (relatively fast, but you can pass them in if you'd like)
    N = length(ψ),                      # SIZE OF STATEVECTOR
    n = length(device),                 # NUMBER OF QUBITS
    m = round(Int, N^(1/n)),            # NUMBER OF LEVELS ON EACH QUBIT
    T = length(pulses[1]),              # TOTAL DURATION OF EVOLUTION

    # CALCULATED VALUES (pre-calculate these and pass them in to speed up optimizations)
    #= The only values required for evolution are: `ΛD`, `UD`, `a_`.
        All others are used only to calculate those.
        Therefore, if you are providing the required values,
            you may set the others to `nothing`.
        Do *NOT* just omit them, though, since doing so invokes expensive calculations.
    =#
    HD= Devices.static_hamiltonian(device, m),  # DEVICE HAMILTONIAN
    ΛU= eigen(HD),                      # `Eigen` STRUCTURE
    ΛD= ΛU.values,                      # DEVICE HAMILTONIAN EIGENVALUES
    UD= ΛU.vectors,                     # DEVICE HAMILTONIAN EIGENVECTORS
    a_= Utils.algebra(n, m, basis=UD),  # LIST OF ROTATED ANNIHILATION OPERATORS
)
    #= NOTE: Personally I feel like we'll want to start in qubit basis someday.
        If that day comes, here is the code...

    # ROTATE INTO DEVICE BASIS
    ψ .= UD' * ψ
                                                =#

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

    # SOLVE THE SYSTEM OF DIFFERENTIAL EQUATIONS
    #= NOTE:
        This method autoselects an algorithm based on I have no idea what,
            meaning I have no idea what the time complexity or accuracy are likely to be,
            or how I should expect them to scale with increasing system size.
        But, it works *pretty* well for the single-qubit case,
            so I'm willing to treat it as a sort of black-box standard.
    =#
    schrodinger = ODEProblem(interaction!, ψ, (0.0, T))
    solution = solve(schrodinger, save_everystep=false)

    # WRITE FINAL SOLUTION TO THE GIVEN STATEVECTOR
    ψ .= solution.u[end]

    # RE-NORMALIZE THIS STATE
    ψ .= ψ / √abs(ψ'*ψ)

    #= NOTE: Personally I feel like we'll want to start in qubit basis someday.
        If that day comes, here is the code...

    # ROTATE *OUT* OF DEVICE BASIS
    ψ .= UD * ψ
                                                =#
end


function evolve!(
    ψ::AbstractVector{<:Number},
    pulses::AbstractVector{<:Pulses.PulseTemplate},
    device::Devices.Device,
    ::Type{Direct};
    numsteps::Integer = 2000,

    # INFERRED VALUES (relatively fast, but you can pass them in if you'd like)
    N = length(ψ),                      # SIZE OF STATEVECTOR
    n = length(device),                 # NUMBER OF QUBITS
    m = round(Int, N^(1/n)),            # NUMBER OF LEVELS ON EACH QUBIT
    T = length(pulses[1]),              # TOTAL DURATION OF EVOLUTION
    t_= range(0,T,numsteps),            # TIME GRID
    Δt= numsteps > 1 ? t_[2]-t_[1] : T, # DURATION OF EACH TIME STEP

    # CALCULATED VALUES (pre-calculate these and pass them in to speed up optimizations)
    #= The only values required for evolution are: `ΛD`, `UD`, `a_`.
        All others are used only to calculate those.
        Therefore, if you are providing the required values,
            you may set the others to `nothing`.
        Do *NOT* just omit them, though, since doing so invokes expensive calculations.
    =#
    HD= Devices.static_hamiltonian(device, m),  # DEVICE HAMILTONIAN
    ΛU= eigen(HD),                      # `Eigen` STRUCTURE
    ΛD= ΛU.values,                      # DEVICE HAMILTONIAN EIGENVALUES
    UD= ΛU.vectors,                     # DEVICE HAMILTONIAN EIGENVECTORS
    a_= Utils.algebra(n, m, basis=UD),  # LIST OF ROTATED ANNIHILATION OPERATORS
)
    ######################################################################################
    #                                 TIME EVOLUTION

    #= NOTE: Personally I feel like we'll want to start in qubit basis someday.
        If that day comes, here is the code...

    # ROTATE INTO DEVICE BASIS
    ψ .= UD' * ψ
                                                =#

    for i ∈ 1:numsteps
        # CONSTRUCT CONTROL HAMILTONIAN (IN DEVICE BASIS)
        HC = zeros(N,N)
        for q ∈ 1:n
            Ω = Pulses.amplitude(pulses[q], t_[i])
            ν = Pulses.frequency(pulses[q], t_[i])
            z = Ω * exp(im*ν*t_[i])
            HC += z*a_[q] + z'*a_[q]'
        end

        # CONJUGATE WITH ACTION OF (DIAGONALIZED) DEVICE HAMILTONIAN
        expD = Diagonal(exp.((im*t_[i]) * ΛD))  # DEVICE ACTION
        HIC = expD * HC * expD'     # INTERACTION-PICTURE CONTROL HAMILTONIAN

        # TODO: pre-allocate HC, expD, and HIC

        # APPLY ACTION OF THE INTERACTION-PICTURE CONTROL HAMILTONIAN
        ψ .= exp( (-im*Δt) * HIC) * ψ
            # TODO: Lanczos is a carbon copy of this method, except this one line.
    end

    #= NOTE: Personally I feel like we'll want to start in qubit basis someday.
        If that day comes, here is the code...

    # ROTATE *OUT* OF DEVICE BASIS
    ψ .= UD * ψ
                                                =#

    ######################################################################################
end


function evolve!(
    ψ::AbstractVector{<:Number},
    pulses::AbstractVector{<:Pulses.PulseTemplate},
    device::Devices.Device,
    ::Type{Rotate};
    numsteps::Integer = 2000,
    qubitapplymode::QubitApplyMode = Kronec(),

    # INFERRED VALUES (relatively fast, but you can pass them in if you'd like)
    N = length(ψ),                      # SIZE OF STATEVECTOR
    n = length(device),                 # NUMBER OF QUBITS
    m = round(Int, N^(1/n)),            # NUMBER OF LEVELS ON EACH QUBIT
    T = length(pulses[1]),              # TOTAL DURATION OF EVOLUTION
    t_= range(0,T,numsteps),            # TIME GRID
    Δt= numsteps > 1 ? t_[2]-t_[1] : T, # DURATION OF EACH TIME STEP

    # CALCULATED VALUES (pre-calculate these and pass them in to speed up optimizations)
    #= The only values required for evolution are: `ΛD`, `UD`, `V`, `a`
        All others are used only to calculate those.
        Therefore, if you are providing the required values,
            you may set the others to `nothing`.
        Do *NOT* just omit them, though, since doing so invokes expensive calculations.
    =#
    HD= Devices.static_hamiltonian(device, m),  # DEVICE HAMILTONIAN
    ΛU= eigen(HD),                      # `Eigen` STRUCTURE
    ΛD= ΛU.values,                      # DEVICE HAMILTONIAN EIGENVALUES
    UD= ΛU.vectors,                     # DEVICE HAMILTONIAN EIGENVECTORS
    V = UD * Diagonal(exp.((-im*Δt) * ΛD)) * UD',   # REPEATED DEVICE ACTION
    a = Utils.a_matrix(m),              # SINGLE-QUBIT ANNIHILATION OPERATOR

    # PRE-ALLOCATIONS (for those that want every last drop of efficiency...)
    O_ = [Matrix{ComplexF64}(undef, m,m) for q ∈ 1:n]   # HOLDS QUBIT-WISE DRIVE OPERATORS
)
    ######################################################################################
    #                                 TIME EVOLUTION

    #= NOTE: The very first step is, mathematically, exp(-𝒊 HD t_[1]),
        but since t_[1]=0, this is an identity operation and we can skip it. =#

    # ROTATE *OUT* OF DEVICE BASIS
    ψ .= UD * ψ
    #= NOTE: Personally I feel like we'll want to start in qubit basis someday.
        If that day comes, delete this first rotation... =#

    # APPLY FIRST PULSE DRIVES  (treated separately to give `V` proper "join" behavior)
    _preparequbitdrives!(pulses, m, t_[1], Δt; n=n, a=a, O_=O_)
    _applyqubitoperators!(ψ, O_, qubitapplymode; N=N, n=n, m=m)

    for i ∈ 2:numsteps
        # CONNECT EACH TIME STEP WITH THE DEVICE ACTION
        ψ .= V * ψ

        # APPLY PULSE DRIVES
        _preparequbitdrives!(pulses, m, t_[i], Δt; n=n, a=a, O_=O_)
        _applyqubitoperators!(ψ, O_, qubitapplymode; N=N, n=n, m=m)
    end

    # LAST STEP: exp(𝒊 HD t[numsteps])), ie. exp(-𝒊 HD T)
    ψ .= UD' * ψ                        # ROTATE INTO DEVICE BASIS
    ψ .*= exp.( (im*T) * ΛD)            # ROTATE PHASES FOR ONE LAST TIME EVOLUTION

    #= NOTE: Personally I feel like we'll want to start in qubit basis someday.
        If that day comes, here is the code...
    ψ .= UD  * ψ                        # ROTATE *OUT* OF DEVICE BASIS
                                                =#

    ######################################################################################
end


"""
    _preparequbitdrives!(
        pulses::AbstractVector{<:Pulses.PulseTemplate}, # PULSE TEMPLATES FOR EACH QUBIT
        m::Integer,                                     # NUMBER OF LEVELS ON EACH QUBIT
        t::Number,                                      # TIME POINT
        Δt::Number;                                     # TIME TO THE NEXT TIME POINT

        # INFERRED VALUES (relatively fast, but you can pass them in if you'd like)
        n = length(pulses),                             # NUMBER OF QUBITS

        # CALCULATED VALUES (pre-calculate these and pass them in to speed up optimizations)
        a = Utils.a_matrix(m),                          # SINGLE-QUBIT ANNIHILATION OPERATOR

        # PRE-ALLOCATIONS (for those that want every last drop of efficiency...)
        O_ = [Matrix{ComplexF64}(undef, m,m) for q ∈ 1:n]   # HOLDS QUBIT-WISE DRIVE OPERATORS
    )


Prepare a vector of qubit operations representing the instantaneous action of a pulse.

Say a pulse has amplitude Ω and frequency ν, and define z = Ω exp(𝒊 ν t)
    We may model the action of the pulse on a resonant system at time t
        with a "Control" Hamiltonian H = z a + z' a',
        and the evolution over short time Δt as exp(-𝒊 Δt H).

"""
function _preparequbitdrives!(
    pulses::AbstractVector{<:Pulses.PulseTemplate}, # PULSE TEMPLATES FOR EACH QUBIT
    m::Integer,                                     # NUMBER OF LEVELS ON EACH QUBIT
    t::Number,                                      # TIME POINT
    Δt::Number;                                     # TIME TO THE NEXT TIME POINT

    # INFERRED VALUES (relatively fast, but you can pass them in if you'd like)
    n = length(pulses),                             # NUMBER OF QUBITS

    # CALCULATED VALUES (pre-calculate these and pass them in to speed up optimizations)
    a = Utils.a_matrix(m),                          # SINGLE-QUBIT ANNIHILATION OPERATOR

    # PRE-ALLOCATIONS (for those that want every last drop of efficiency...)
    O_ = [Matrix{ComplexF64}(undef, m,m) for q ∈ 1:n]   # HOLDS QUBIT-WISE DRIVE OPERATORS
)
    ######################################################################################
    for q ∈ 1:n
        # EXTRACT TIME-DEPENDENT COEFFICIENTS
        Ω = Pulses.amplitude(pulses[q], t)
        ν = Pulses.frequency(pulses[q], t)
        z = Ω * exp(im*ν*t)

        # CONSTRUCT AND EXPONENTIATE MATRIX
        O_[q] .= exp((-im*Δt) * (z*a + z'*a'))
    end
    ######################################################################################
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
This will simply combine all qubit operators into a full NxN matrix and left-multiply ψ.

"""
function _applyqubitoperators!(
    ψ::AbstractVector{<:Number},
    O_::AbstractVector{<:AbstractMatrix{<:Number}},
    mode::QubitApplyMode;

    # INFERRED VALUES (relatively fast, but you can pass them in if you'd like)
    #= None of these are needed for mode::Kronec, you may pass them in as `nothing` =#
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


function evolve!(
    ψ::AbstractVector{<:Number},
    pulses::AbstractVector{<:Pulses.PulseTemplate},
    device::Devices.Device,
    ::Type{Prediag};
    numsteps::Integer = 2000,
    suzukiorder = 2,
    qubitapplymode::QubitApplyMode = Kronec(),

    # INFERRED VALUES (relatively fast, but you can pass them in if you'd like)
    N = length(ψ),                      # SIZE OF STATEVECTOR
    n = length(device),                 # NUMBER OF QUBITS
    m = round(Int, N^(1/n)),            # NUMBER OF LEVELS ON EACH QUBIT
    T = length(pulses[1]),              # TOTAL DURATION OF EVOLUTION
    t_= range(0,T,numsteps),            # TIME GRID
    Δt= numsteps > 1 ? t_[2]-t_[1] : T, # DURATION OF EACH TIME STEP

    # CALCULATED VALUES (pre-calculate these and pass them in to speed up optimizations)
    #= The only values required for evolution are:
            `ΛD`, `UD`, `L`, `Λ`, `UQP`, *`UPQ`, `in_basis`, `outbasis`
                *Even this one is not required for `suzukiorder=1`.
        All others are used only to calculate those.
        Therefore, if you are providing the required values,
            you may set the others to `nothing`.
        Do *NOT* just omit them, though, since doing so invokes expensive calculations.
    =#
    HD= Devices.static_hamiltonian(device, m),  # DEVICE HAMILTONIAN
    ΛU= eigen(HD),                      # `Eigen` STRUCTURE
    ΛD= ΛU.values,                      # DEVICE HAMILTONIAN EIGENVALUES
    UD= ΛU.vectors,                     # DEVICE HAMILTONIAN EIGENVECTORS
    V = UD * Diagonal(exp.((-im*Δt) * ΛD)) * UD',   # REPEATED DEVICE ACTION

    a = Utils.a_matrix(m),              # SINGLE-QUBIT ANNIHILATION OPERATOR
    Q =      (a + a'),                  # CANONICAL COORDINATE OPERATOR
    P = im * (a - a'),                  # CANONICAL   MOMENTUM OPERATOR
    ΛUQ = eigen(Q),
    ΛUP = eigen(P),
    Λ = ΛUQ.values,                     # EIGENVALUES OF Q OPERATOR (OR P)
    UQ= ΛUQ.vectors,
    UP= ΛUP.vectors,
    UQP = UQ' * UP,                     # ROTATION MATRIX FROM P -> Q BASIS
    UPQ = UP' * UQ,                     # ROTATION MATRIX FROM Q -> P BASIS

    in_basis = Utils.kron_concat(suzukiorder==1 ? UP : UQ, n),  # STARTING BASIS FOR DRIVE
    outbasis = Utils.kron_concat(UQ, n),                        #   ENDING BASIS FOR DRIVE
    L = in_basis' * V * outbasis,       # LIGAND OPERATION

    # PRE-ALLOCATIONS (for those that want every last drop of efficiency...)
    O_ = [Matrix{ComplexF64}(undef, m,m) for q ∈ 1:n]   # HOLDS QUBIT-WISE DRIVE OPERATORS
)
    ######################################################################################
    #                                 TIME EVOLUTION

    #= NOTE: The very first step is, mathematically, exp(-𝒊 HD t_[1]),
        but since t_[1]=0, this is an identity operation and we can skip it. =#

    # ROTATE *OUT* OF DEVICE BASIS
    ψ .= UD * ψ
    #= NOTE: Personally I feel like we'll want to start in qubit basis someday.
        If that day comes, delete this first rotation... =#

    # ROTATE INTO `in_basis`
    ψ .= in_basis' * ψ

    # APPLY FIRST PULSE DRIVES  (treated separately to give `V` proper "join" behavior)
    _preparequbitdrives_productformula(pulses, m, t_[1], Δt; suzukiorder=suzukiorder,
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

    # ROTATE *OUT* OF `outbasis`
    ψ .= outbasis * ψ

    # LAST STEP: exp(𝒊 HD t[numsteps])), ie. exp(-𝒊 HD T)
    ψ .= UD' * ψ                        # ROTATE INTO DEVICE BASIS
    ψ .*= exp.( (im*T) * ΛD)            # ROTATE PHASES FOR ONE LAST TIME EVOLUTION

    #= NOTE: Personally I feel like we'll want to start in qubit basis someday.
        If that day comes, here is the code...
    ψ .= UD  * ψ                        # ROTATE *OUT* OF DEVICE BASIS
                                                =#

    ######################################################################################
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

Say a pulse has amplitude Ω and frequency ν, and define z = Ω exp(𝒊 ν t)
    We may model the action of the pulse on a resonant system at time t
        with a "Control" Hamiltonian H = z a + z' a',
        and the evolution over short time Δt as exp(-𝒊 Δt H).

This variant rewrites the drive Hamiltonian (H = z a + z' a') => x Q + y P,
    to rewrite the evolution operator exp(-𝒊 Δt H) ≈ exp(-𝒊 Δt x Q) exp(-𝒊 Δt x P)
        or a related product formula, selected with `suzukitrotter`.

"""
function _preparequbitdrives_productformula(
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

    # # CALCULATED VALUES (pre-calculate these and pass them in to speed up optimizations)
    # a = Utils.a_matrix(m),                          # SINGLE-QUBIT ANNIHILATION OPERATOR
    # Q =      (a + a'),                              # CANONICAL COORDINATE OPERATOR
    # P = im * (a - a'),                              # CANONICAL   MOMENTUM OPERATOR
    # ΛUQ = eigen(Q),
    # ΛUP = eigen(P),
    # Λ = ΛUQ.values,
    # UQ= ΛUQ.vectors,
    # UP= ΛUP.vectors,
    # UQP = UQ * UP',
    # UPQ = UP * UQ',

    # PRE-ALLOCATIONS (for those that want every last drop of efficiency...)
    O_ = [Matrix{ComplexF64}(undef, m,m) for q ∈ 1:n]   # HOLDS QUBIT-WISE DRIVE OPERATORS
)
    ######################################################################################
    for q ∈ 1:n
        # EXTRACT TIME-DEPENDENT COEFFICIENTS
        Ω = Pulses.amplitude(pulses[q], t)
        ν = Pulses.frequency(pulses[q], t)
        z = Ω * exp(im*ν*t)
        x, y = real(z), imag(z)

        # EVOLVE QUBIT IN TIME, AND EXTEND FULL-QUBIT OPERATOR
        if     suzukiorder == 1
            expQ = Diagonal(exp.((-im*Δt*real(z)) * Λ))
            expP = Diagonal(exp.((-im*Δt*imag(z)) * Λ))

            O_[q] .= expQ * UQP * expP
        elseif suzukiorder == 2
            expQ = Diagonal(exp.((-im*Δt*real(z)/2) * Λ))
            expP = Diagonal(exp.((-im*Δt*imag(z)  ) * Λ))

            O_[q] .= expQ * UQP * expP * UPQ * expQ
        else
            error("Only `suzukiorder`s 1 and 2 are supported.")
        end
    end

    # TODO: Pre-allocate expQ, expP

    ######################################################################################
    return O_
end

end
