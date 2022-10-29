#= Check direct evolution Julia implementation against the one in ctrlq. =#

include("../src/utils.jl")
include("../src/device.jl")
include("../src/pulse.jl")
include("../src/evolve.jl")

using LinearAlgebra
import ..Devices
import ..Pulses
import ..Evolutions

##########################################################################################
#                               SCALAR CONTROL PARAMETERS
m = 3                   # NUMBER OF LEVELS ON EACHQUBIT
n = 4                   # NUMBER OF QUBITS
T = 10.0                # DURATION (ns)
numsteps = 2000         # NUMBER OF TIME DIVISIONS

##########################################################################################
#                               CONTROL OBJECTS

# CONSTRUCT INITIAL STATE
kets = [1,1,1,1,1,1][1:n]   # SPECIFY INITIAL BASIS VECTOR (gets truncated to n qubits)

_init = [1]
for q ∈ 1:n
    ψq = zeros(m);  ψq[kets[q]+1] = 1
    global _init = kron(_init, ψq)
end
ψI::Vector{Complex} = _init


# DESIGN PULSE PATTERNS
pulses = [
    Pulses.BasicSquarePulse(T, 29.0, [0.5, 0.4], [2.3]),
    Pulses.BasicSquarePulse(T, 31.0, [0.3, 0.2], [5.6]),
    Pulses.BasicSquarePulse(T, 29.2, [0.3, 0.1], [6.1]),
    Pulses.BasicSquarePulse(T, 29.8, [0.5, 0.6], [1.0]),
    Pulses.BasicSquarePulse(T, 30.6, [0.4, 0.2], [3.4]),
    Pulses.BasicSquarePulse(T, 31.1, [0.5, 0.1], [5.6]),
][1:n]  # (truncate to n qubits)

# DESIGN DEVICE
device = Devices.Transmon(
    2π .* [
        4.808049015463495,
        4.833254817254613,
        4.940051121317842,
        4.795960998582043,
        5.327178905327263,
        5.076070923892718,
    ],
    2π .* [
        0.018312874435769682,
        0.021312874435769682,
        0.019312874435769682,
        0.020312874435769682,
        0.310177361313422900,
        0.291617038572545600,
    ],
    Dict(
        Devices.QubitCouple(1,2) => 0.1,
        Devices.QubitCouple(1,4) => 0.2,
        Devices.QubitCouple(1,6) => 0.3,
        Devices.QubitCouple(2,3) => 0.4,
        Devices.QubitCouple(2,5) => 0.5,
        Devices.QubitCouple(3,4) => 0.6,
        Devices.QubitCouple(3,6) => 0.7,
        Devices.QubitCouple(4,5) => 0.8,
        Devices.QubitCouple(5,6) => 0.9,
    ),
    n,  # (truncate to n qubits)
)


N = m^n                                     # DIMENSION OF HILBERT SPACE
HD = Devices.static_hamiltonian(device, m)  # STATIC HAMILTONIAN

##########################################################################################
#                                     ctrlq SETUP
using PyCall

# CONSTRUCT PULSE OBJECT
ctrlq_solve = pyimport("ctrlq.lib.solve")
pobjc = ctrlq_solve.pulsec(
    [pulse.amplitudes for pulse ∈ pulses],
    [pulse.steptimes  for pulse ∈ pulses],
    [pulse.frequency  for pulse ∈ pulses],
    T, n, length(pulses[1].amplitudes),
)

# CONSTUCT TRANSMON OBJECT
cvqe = pyimport("ctrlq.cvqe")
myham = cvqe.transmon(
    nqubit=n,
    nstate=m,
    mham=zeros(1,1),
    # istate=kets,  # ...wish this worked..!
    Hstatic=HD,
)
myham.initialize_psi(kets)

#######      FUNCTION HANDLER
cvqe_evolve = pyimport("ctrlq.cvqe.evolve")
function ctrlqevolve()
    return cvqe_evolve.evolve(
        ψI,
        pobjc,
        myham,
        solver="trotter",
        nstep=numsteps,
    )[:,1]
end




##########################################################################################
#                               PRE-CALCULATE INITIAL SETUP


#= DRESSED BASIS

My code uses the "device basis" specified by the eigenvectors returned by eigen(...).
These eigenvectors are ordered in ascending order of the corresponding eigenvalues,
    and I think their phase is normalized such that the first row is real and positive.

Meanwhile, ctrlq orders eigenvectors such that the ith component of the ith eigenvector
        is the largest in that vector,
    and the phase is is normalized such that that component is real and positive.
This ordering makes the rotation matrix look more similar to an identity operation,
    which I'm sure has nice properties.
But, it strikes me as unnecessarily complicated, so I'm not using it in my code.

Alas, that means that any direct comparison between statevectors in the two codes
    must introduce this generalized permutation matrix one way or another.
In this script, I'll simply pre-calculate the "dressed" basis ctrlq uses,
    and force my code to use it instead of the default.

=#

# DIAGONALIZE THE DEVICE HAMILTONIAN
ΛD, UD = eigen(HD)                          # DIAGONALIZED

# IMPOSE PERMUTATION
σ = Vector{Int}(undef,N)
for i in 1:N
    perm = sortperm(abs.(UD[i,:]), rev=true)# STABLE SORT BY "ith" COMPONENT
    perm_ix = 1                             # CAREFULLY HANDLE TIES
    while perm[perm_ix] ∈ σ[1:i-1]
        perm_ix += 1
    end
    σ[i] = perm[perm_ix]
end
ΛD .= ΛD[  σ]
UD .= UD[:,σ]

# IMPOSE PHASE
for i in 1:N
    if UD[i,i] < 0; UD[:,i] *= -1;  end
end

# DOUBLE-CHECK FOR CORRECTNESS
@assert HD ≈ (UD * Diagonal(ΛD) * UD')

#= TIME STEP

ctrlq creates the range `t=linspace(0, T, N=numsteps)`,
    then it calculates the spacing as `Δt=T/numsteps`.
Actually, this is the spacing for `linspace(0,T,N=numsteps,endpoint=False)`,
    or better yet `linspace(0,T,N=numsteps+1)[1:]`,
    which is actually the range I think maybe we *should* use.
But for now, the correct spacing is `T/(numsteps-1)`.

That said, since we want to confirm our Julia code matches ctrlq,
    we'll manually force the erroneous Δt.
=#
t = range(0,T,numsteps)
Δt = T / numsteps
# Δt = t[2]-t[1]      # UNCOMMENT TO SEE HOW MUCH CORRECT Δt CHANGES RESULTS

a_= Utils.algebra(n,m,basis=UD)             # ANNIHILATION OPERATORS ON EACH QUBIT

#######      FUNCTION HANDLER
function juliaevolve()
    ψ = copy(ψI)
    Evolutions.evolve!(
        ψ, pulses, device, Evolutions.Direct;
        numsteps=numsteps,
        N=N, n=n, m=m,
        T=T, t_=t, Δt=Δt,
        HD=nothing, ΛU=nothing,
        ΛD=ΛD, UD=UD, a_=a_,
    )
    return ψ
end




##########################################################################################
#                                    COMPARE RUNS

fidelity(ψ,φ) = 1 - abs2(ψ'*φ)

ψJ = juliaevolve()
println("1-|⟨ψJ|ψJ⟩|² = $(fidelity(ψJ,ψJ))")

ψP = ctrlqevolve()
println("1-|⟨ψP|ψP⟩|² = $(fidelity(ψP,ψP))")

println("1-|⟨ψJ|ψP⟩|² = $(fidelity(ψJ,ψP))")
#=

Works perfect for my n=2 test case.

Not so good for my n=4 test case.
    Why? Well, unfortunately, the ctrlq Hamiltonian doesn't seem to match mine after dressing.
    Is my method of permuting not working?!
    No, I think it's working just fine. I think that ctrlq's method is broken.

    Presumptuous, no? Well, I'm pretty sure that myham.dsham is *supposed* to be diagonal. It isn't.
    And even if I'm wrong about the idea that dsham should be the device hamiltonian in the diagonal dressed basis, *surely* it should at least be the device hamiltonian, no? So it will have the same eigenvalues as HD. Ahm, it doesn't.
=#
eigsP = eigen(myham.dsham.toarray()).values
eigs  = eigen(HD).values
@assert eigsP ≈ eigs





# using BenchmarkTools
# println("Julia benchmarking:")
# @btime juliaevolve()
#
# println("ctrlq benchmarking:")
# @btime ctrlqevolve()
