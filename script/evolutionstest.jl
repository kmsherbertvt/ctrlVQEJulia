#= Validate that each evolution method gives roughly the same results. =#

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
ψI::Vector{ComplexF64} = _init


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


##########################################################################################
#                           PERFORM THE ACTUAL EVOLUTIONS

fidelity(ψ,φ) = 1 - abs2(ψ'*φ)

# DIRECT EXPONENTIATION (Checked against ctrlq, this is the "standard".)
ψD = Evolutions.evolve(ψI, pulses, device, Evolutions.Direct; numsteps=numsteps)
println("1-|⟨ψD|ψI⟩|²: $(fidelity(ψD, ψI))")        # CHECK HOW FAR WE EVOLVED FROM ψI

# NUMERICAL INTEGRATION
ψS = Evolutions.evolve(ψI, pulses, device, Evolutions.ODE; numsteps=numsteps)
println("1-|⟨ψD|ψS⟩|²: $(fidelity(ψD, ψS))")

# ROTATE BETWEEN STATIC AND DRIVE BASES
ψRk = Evolutions.evolve(
    ψI, pulses, device, Evolutions.Rotate; numsteps=numsteps,
    qubitapplymode=Evolutions.Kronec(), # APPLY QUBIT ROTATIONS WITH KRONECKER PRODUCT
)
println("1-|⟨ψD|ψRk⟩|²: $(fidelity(ψD, ψRk))")

ψRt = Evolutions.evolve(
    ψI, pulses, device, Evolutions.Rotate; numsteps=numsteps,
    qubitapplymode=Evolutions.Tensor(), # APPLY QUBIT ROTATIONS WITH TENSOR ALGEBRA
)
println("1-|⟨ψD|ψRt⟩|²: $(fidelity(ψD, ψRt))")

# FACTOR DRIVE BASIS SO ALL ROTATIONS ARE TIME-INDEPENDENT
ψ1k = Evolutions.evolve(
    ψI, pulses, device, Evolutions.Prediag; numsteps=numsteps,
    suzukiorder=1,                      # FACTOR WITH SIMPLEST POSSIBLE PRODUCT FORMULA
    qubitapplymode=Evolutions.Kronec(), # APPLY QUBIT ROTATIONS WITH KRONECKER PRODUCT
)
println("1-|⟨ψD|ψ1k⟩|²: $(fidelity(ψD, ψ1k))")

ψ1t = Evolutions.evolve(
    ψI, pulses, device, Evolutions.Prediag; numsteps=numsteps,
    suzukiorder=1,                      # FACTOR WITH SIMPLEST POSSIBLE PRODUCT FORMULA
    qubitapplymode=Evolutions.Tensor(), # APPLY QUBIT ROTATIONS WITH TENSOR ALGEBRA
)
println("1-|⟨ψD|ψ1t⟩|²: $(fidelity(ψD, ψ1t))")

ψ2k = Evolutions.evolve(
    ψI, pulses, device, Evolutions.Prediag; numsteps=numsteps,
    suzukiorder=2,                      # FACTOR WITH SYMMETRIC PRODUCT FORMULA
    qubitapplymode=Evolutions.Kronec(), # APPLY QUBIT ROTATIONS WITH KRONECKER PRODUCT
)
println("1-|⟨ψD|ψ2k⟩|²: $(fidelity(ψD, ψ2k))")

ψ2t = Evolutions.evolve(
    ψI, pulses, device, Evolutions.Prediag; numsteps=numsteps,
    suzukiorder=2,                      # FACTOR WITH SYMMETRIC PRODUCT FORMULA
    qubitapplymode=Evolutions.Tensor(), # APPLY QUBIT ROTATIONS WITH TENSOR ALGEBRA
)
println("1-|⟨ψD|ψ2t⟩|²: $(fidelity(ψD, ψ2t))")
