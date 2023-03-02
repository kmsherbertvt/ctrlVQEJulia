#= Shotgun search over initial parameter space, recording what energy they optimize to. =#


##########################################################################################

# WHICH SEEDS TO RUN
sd_ = 1:100

# PULSE HYPER-PARAMETERS
T = 100.0       # TOTAL PULSE DURATION
W = 200         # NUMBER OF PULSE WINDOWS

# THE NUMBER OF MODES CONSIDERED FOR EACH TRANSMON
m = 2

# NUMBER OF TIME STEPS
r = 1000

# PULSE CONSTRAINTS
ΔΩ = 2π * 0.02      # RANGE OF PULSE AMPLITUDE
ΔΔ = 2π * 1.0       # RANGE OF DE-TUNING

##########################################################################################

# `using` STATEMENTS
using Random, LinearAlgebra         # STANDARD LIBRARIES
using NPZ, LBFGSB                   # NEED TO BE INSTALLED

# `include` STATEMENTS
include("../src/utils.jl")
include("../src/device.jl")
include("../src/pulse.jl")
include("../src/evolve.jl")
include("../src/gradient.jl")
include("../chemistry/molecule.jl")

# `import` STATEMENTS
import ..Utils, ..Devices, ..Pulses, ..Evolutions, ..Gradients      # CORE MODULES
import ..Molecules                                                  # CHEMISTRY "PLUG-IN"

##########################################################################################

# LOAD THE HAMILTONIAN FROM AN `.npy` FILE
Hname = "h6ham15"
H = npzread("matrix/$Hname.npy")
n = round(Int, log2(size(H,1)))

# PHYSICAL CHARACTERISTICS OF OUR DEVICE
# device = Devices.selectqubits(1:n, Devices.Transmon(
#     2π*[4.8080, 4.8333, 4.9400, 4.7960],        # QUBIT RESONANCE FREQUENCIES
#     2π*[0.3102, 0.2916, 0.3302, 0.2616],        # QUBIT ANHARMONICITIES
#     Dict{Devices.QubitCouple,Float64}(          # QUBIT COUPLING CONSTANTS
#         Devices.QubitCouple(1,2) => 2π*.01831,
#         Devices.QubitCouple(2,3) => 2π*.02131,
#         Devices.QubitCouple(3,4) => 2π*.01931,
#         Devices.QubitCouple(4,1) => 2π*.02031,
#     )
# ))
# devicename = "npj"

ω = 2π*collect(4.8 .+ (0.02 * (1:n)))
δ = 2π*0.3 * ones(n)
G = Dict{Devices.QubitCouple,Float64}()
for p in 1:n
    q = (p == n) ? 1 : p + 1
    G[Devices.QubitCouple(p,q)] = 2π*0.02
end
device = Devices.Transmon(ω, δ, G)
devicename = "sys"

# FILENAME TO WRITE TO. BE SURE TO CHANGE THIS IF THE OBSERVABLE OR DEVICE CHANGE
filename = "dat/squarefixν.scan.$Hname.$devicename"

# WRITE HEADER IF THE FILE IS UNUSED
if !isfile("$filename.csv")
    open("$filename.csv", "w") do io
        println(io, join([
            "T", "W", "ΔΩ", "ΔΔ", "m", "r", "sd",
            "E0", "rt", "Ex", "lk", "En", "εE", "fC"
        ], "\t"))

    end
end

# FETCH THE BEST TRIAL FOUND SO FAR
npzfilename = "$filename.opt.T($T).npz"
min_En = isfile(npzfilename) ? npzread(npzfilename)["En"] : Inf

##########################################################################################

# TOTAL SIZE OF HILBERT SPACE BEING SIMULATED
N = m^n

# SINGLE-QUBIT BOSONIC ANNIHILATION OPERATOR (truncated to m modes)
a = Utils.a_matrix(m)

# EXTEND OBSERVABLE ONTO SIMULATION SPACE
Π = Utils.projector(n, 2, m)    # PROJECTOR FROM SIMULATION SPACE ONTO PROBLEM SPACE
O = Hermitian(Π'*H*Π)           # MOLECULAR HAMILTONIAN, EMBEDDED INTO SIMULATION SPACE

# PREPARE THE HARTREE-FOCK STATE
ψ0 = zeros(ComplexF64, N)                               # INITIALIZE ⟨z|ψ0⟩ = 0
ψ0[argmin(real.(diag(O)))] = one(ComplexF64)            #    ASSIGN ⟨HF|ψ0⟩ = 1

# SET OPTIONS FOR DIFFERENT CALCULATION MODES
qubitapplymode = Evolutions.Kronec()    # HOW TO APPLY SINGLE-QUBIT OPERATORS?
iobasis = Evolutions.QubitBasis()       # HOW TO INTERPRET THE INPUT/OUTPUT STATEVECTOR?

# DISCRETIZE TIME
t_= range(0,T,r+1)                 # TIME GRID
τ = T / r                          # DURATION OF EACH TIME STEP

# HEAVY CALCULATIONS FOR WORKING IN THE DEVICE BASIS (aka "dressed basis")
HD = Devices.static_hamiltonian(device, m)  # DEVICE HAMILTONIAN
ΛD, UD = Utils.dressedbasis(HD)             # DEVICE EIGENVALUES AND EIGENVECTORS
L = UD* Diagonal(exp.((-im*τ) * ΛD)) *UD'   # REPEATED DEVICE ACTION FOR EACH TIME STEP

# PRE-ALLOCATIONS
ψ  = Vector{ComplexF64}(undef, m^n)         # FOR STORING THE WAVEFUNCTION
∇Ω = Matrix{Float64}(undef, r+1, n)         # FOR STORING THE DERIVATIVES

_N1 = Vector{ComplexF64}(undef, N)          # FOR MATRIX-VECTOR MULTIPLICATION
_m2_ = [Matrix{ComplexF64}(undef, m,m) for q ∈ 1:n]     # QUBIT-WISE DRIVE OPERATORS
_K_ = [Matrix{ComplexF64}(undef, m^q, m^q) for q ∈ 1:n] # FOR KRON'ing OPERATORS

# WRAP ALL THESE UP INTO A NICE, TIDY KEYWORD DICTIONARY
evolve_kwargs = Dict(
    :iobasis => iobasis, :qubitapplymode => qubitapplymode, :numsteps => r,
    :N => N, :n => n, :m => m, :T => T, :t_=> t_, :Δt=> τ,
    :ΛD => ΛD, :UD => UD, :V  => L, :a  => a,
    :tmpV  => _N1, :tmpM_ => _m2_, :tmpK_ => _K_,
)

gradient_kwargs = Dict(
    :iobasis => iobasis, :qubitapplymode => qubitapplymode, :r => r,
    :N => N, :n => n, :m => m, :T => T, :t_=> t_, :τ => τ,
    :ΛD => ΛD, :UD => UD, :V  => L, :a  => a,
    :tmpV  => _N1, :tmpM_ => _m2_, :tmpK_ => _K_, :∂Ω => ∇Ω,
)

# PULSE WINDOW STEP TIMES
steptimes = repeat( range(0,T,W+1)[2:end-1], 1, n)

# CONFIGURE BOUNDS
lbounds = Matrix{Float64}(undef, W, n); ubounds = Matrix{Float64}(undef, W, n)
# # #
lbounds[:,:] .= -ΔΩ                       # SET MAXIMUM AMPLITUDE
ubounds[:,:] .=  ΔΩ

# PRE-ALLOCATE SPACE FOR PARAMETERS (so we can look at them outside loop)
ν  = Vector{Float64}(undef, n)
x0 = Vector{Float64}(undef, W*n)
x  = Vector{Float64}(undef, W*n)

##########################################################################################

"""
    constructpulse(x::AbstractVector{<:Real})

This function constructs a set of pulse objects for each qubit,
    out of a vector of optimizable parameters.

In this example, each pulse will be a sequence of `W` square pulses.

The vector `x` must specify:
- `W` different amplitudes A[j]
- a single ν, for the frequency of the pulse.

Actually, `x` must specify each of these for all `n` qubits,
    for a total of `n(W+1)` parameters.

"""
function constructpulses(x::AbstractVector{<:Real})
    x = reshape(x, :, n)    # RE-INTERPRET `x` AS `n` COLUMNS

    pulses = Vector{Pulses.BasicSquarePulse}(undef, n)
    for q in 1:n
        A = x[:,q]                # AMPLITUDE ARRAY
        pulses[q] = Pulses.BasicSquarePulse(T, ν[q], A, steptimes[:,q])
    end
    return pulses
end


"""
    pulsegradient_amplitude(t, j, pulse)

Consider the pulse amplitude Ω(t) on a single qubit.
    In our example, Ω(t) is a step function, determined by each A[j] and s[j].
    We're assuming s[j] are held constant in this tutorial,
        but a bit of calculus gives us the partial derivatives ∂Ω(t)/∂A[j].

This function gives ∂Ω(t)/∂A[j].

"""
function pulsegradient_amplitude(t, j, pulse)
    return Utils.interval(t,
        get(pulse.steptimes, j-1, 0),               #  LEFT BOUND (TREAT s[0] as t=0)
        get(pulse.steptimes, j, pulse.duration),    # RIGHT BOUND (TREAT s[W] as t=T)
    )
end

∂Ω∂x = Vector{Float64}(undef, r+1)  # PRE-ALLOCATE A VECTOR TO HOLD THESE GRADIENTS
h = 1e-6                # SPACING FOR A FINITE DIFFERENCE ON THE FREQUENCY PARAMETER


"""
    f(x)

The cost function to be used in optimization.

Returns the energy ⟨ψ|O|ψ⟩, where |ψ⟩ is the quantum state prepared
    by applying the pulse sequences parameterized by `x` to the initial state |ψ0⟩.

Note that this energy is "unnormalized",
    meaning it will actually be implicitly a little higher
    than the one you'd get by ⟨ψ'|H|ψ'⟩, where |ψ'⟩ is |ψ⟩ retricted to 0's and 1's.
That means leakage into modes larger than 0 and 1 will be implicitly penalized
    in an optimization.

The true minimum of O will be exactly equal to the true minimum of H,
    but there's no guarantee an optimization *finds* the true minimum,
    so you must correct for leakage before interpreting the accuracy of your optimization:

    ⟨ψ'|H|ψ'⟩ = ⟨ψ|O|ψ⟩ / ⟨ψ|Π'Π|ψ⟩

"""
function f(x)
    pulses = constructpulses(x)
    ψ .= Evolutions.evolve(ψ0, pulses, device, Evolutions.Rotate; evolve_kwargs...)
    return real(Utils.expectation(O, ψ))
end


"""
    g!(G, x)

The gradient function to be used in optimization.

`G` is just a pre-allocated block of memory to write the gradient vector to.
The `!` in the name is just Julia convention to mark this function as one
    which mutates its primary operand (in this case, `G`).

To calculate the gradients ∂E/∂A[j],
    we'll need to do some rather unrigorous calculus:

    ∂E/∂x[j] = ∂E/∂Ω(t) · ∂Ω(t)/∂x[j]
                        ↑
            This is a vector dot product.

We've implemented ∂Ω(t)/∂x[j] at the start of this section.
The ∂E/∂Ω(t) is what we're calling the gradient signal, or switching function φ(t).

    (Actually my gradient signal has an extra `τ` (aka `dt`) in it,
        and I'm not exactly certain what it's doing there.
    But it seems to work... )

To calculate the gradients ∂E/∂ν...
    ah, well, there is a *frequency* gradient signal I could calculate.
But I haven't implemented that yet.
So what we're going to do instead is run a forward finite difference,
    on just these `n` parameters... ^_^

"""
function g!(G, x)
    G = reshape(G, :, n)    # RE-INTERPRET `G` AS `n` COLUMNS

    # CONSTRUCT THE AMPLITUDE GRADIENT SIGNAL ∇Ω ~= ∂E/∂Ω(t) dt
    pulses = constructpulses(x)
    Gradients.gradientsignal(ψ0, pulses, device, O; gradient_kwargs...)
        #= Because we put a pre-allocated `∇Ω` into `gradient_kwargs`,
            this method writes the gradient signal to that `∇Ω`. =#

    # EVALUATE THE COST FUNCTION AT x, FOR FORWARD FINITE DIFFERENCES
    for (q, pulse) in enumerate(pulses)
        for j in 1:W
            # SET THE GRADIENTS G[j] ~= ∂E/∂A[j]
            ∂Ω∂x .= map(t -> pulsegradient_amplitude(t,j,pulse), t_)
            G[j,q] = ∇Ω[:,q] · ∂Ω∂x
        end
    end

end


##########################################################################################

# GET HARTREE FOCK AND FCI ENERGY
Λ, U = eigen(H)
E_FCI = Λ[1]
ftime = @elapsed E_HF = f(zero(x0))
E_CORR = E_HF - E_FCI

println()
println("*"^50)
println()
println("                FCI energy: $E_FCI")
println("                 HF energy: $E_HF")
println("        Correlation energy: $E_CORR")
println("Time to run single fn evel: ~ $ftime s")
println()

##########################################################################################


function runoptimization(x0)
    # RUN THE OPTIMIZATION
    _, x = lbfgsb(f, g!, x0, lb=vec(lbounds), ub=vec(ubounds),
        m=5,
        factr=1e1,      # UNITS OF MACHINE EPSILON
        pgtol=1e-5,
        iprint=50,
        maxfun=10000,
        maxiter=10000,
    )

    return x
end


# START TAKING DATA
for sd in sd_
    # CONFIGURE INITIAL PARAMETERS
    Random.seed!(sd)
    global ν = ΔΔ * (2*rand(n) .- 1) .+ device.ω
    global x0 = reshape(x0, W, n)
    x0[:,:] .= ΔΩ * (2*rand(W,n) .- 1)
    x0 = vec(x0)

    # INITIAL RUN TO GET A FEEL FOR WHERE WE ARE
    E0 = f(x0)

    println()
    println("*"^50)
    println()
    println("Seed: $sd")
    println("            Initial energy: $E0")
    println(" Distance from FCI, wrt HF: $((E0-E_FCI)/E_CORR)")
    println()

    # RUN OPTIMIZATION
    rt = @elapsed global x .= runoptimization(x0)

    # CONSTRUCT A BUNCH OF USEFUL VALUES FROM OPTIMIZATION RESULTS
    Ex = f(x)                   # CONSTRUCT FINAL ENERGY AND STATEVECTOR
    lk = 1 - real(ψ'*Π'*Π*ψ)    # CALCULATE LEAKAGE
    En = Ex / (1 - lk)          # CALCULATE RE-NORMALIZED ENERGY
    εE = En - E_FCI             # CALCULATE ERROR FROM FCI
    fC = 1 - εE / E_CORR        # CALCULATE FRACTION OF CORRELATION ENERGY RECOVERED

    # WRITE RECORD
    open("$filename.csv", "a") do io
        println(io, join((T, W, ΔΩ, ΔΔ, m, r, sd, E0, rt, Ex, lk, En, εE, fC), "\t"))
    end

    # REPLACE "BEST-CASE" NPZ FILE
    if En < min_En
        npzwrite(npzfilename, Dict(
            "T"   => T,
            "W"   => W,
            "DA"  => ΔΩ,
            "DD"  => ΔΔ,
            "m"   => m,
            "r"   => r,
            "sd"  => sd,
            "E0"  => E0,
            "rt"  => rt,
            "Ex"  => Ex,
            "lk"  => lk,
            "En"  => En,
            "eE"  => εE,
            "fC"  => fC,

            "H"   => H,
            "ω"   => device.ω,
            "δ"   => device.δ,
            "G"   => Devices.getcouplingmatrix(device),
            "n"   => n,
            "y0"  => ψ0,

            "ν"   => ν,
            "x0"  => x0,
            "x"   => x,
        ))
    end
end