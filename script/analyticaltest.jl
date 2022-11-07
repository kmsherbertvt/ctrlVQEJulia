#= Compare each evolution method against the exact solution for a single-body device. =#

include("../src/utils.jl")
include("../src/device.jl")
include("../src/pulse.jl")
include("../src/evolve.jl")
include("../src/analytical.jl")

using LinearAlgebra
import ..Devices
import ..Pulses
import ..Evolutions
import ..Analytical


##########################################################################################
#                                   CONTROL PARAMETERS

ω₀ = 2                  # DEVICE RESONANCE FREQUENCY
δ  = 0                # DEVICE ANHARMONICITY
ν  = 3                # PULSE FREQUENCY
Ω₀ = 1                # PULSE AMPLITUDE
T  = 2.0                # PULSE DURATION (ns)

# DEFINE INITIAL WAVEFUNCTION
# NOTE: Provide three amplitudes here; the third will be ignored for qubit test.
# NOTE: The state is manually normalized later, so provide any three complex numbers here.
#       Actually don't provide [0, 0, ?] since that will make qubit test unnormalizable...
ψI = [
    .3,
    √(1-.3^2)+0im,
    0,  # MANUALLY IGNORED IN QUBIT TEST, NO NEED TO OMIT
]



##########################################################################################
#                               SCALAR CONTROL PARAMETERS

m = 3                   # NUMBER OF LEVELS ON EACHQUBIT
n = 4                   # NUMBER OF QUBITS
T = 10.0                # DURATION (ns)
numsteps = 2000         # NUMBER OF TIME DIVISIONS

##########################################################################################
#                               CONTROL OBJECTS

# DESIGN PULSE PATTERNS
pulses = [Pulses.BasicSquarePulse(T, ν, [Ω₀], Vector{Number}())]

# DESIGN DEVICE
device = Devices.Transmon([ω₀], [δ])

# INITIAL STATES (NORMALIZED)
normalize(ψ) = ψ / √abs(ψ'*ψ)
ψI2 = normalize(ψI[1:2])
ψI3 = normalize(ψI[1:3])

# ANALYTICAL RESULTS
ψA2 = Analytical.onequbitsquarepulse( ψI2, ν, Ω₀, T, ω₀   )
ψA3 = Analytical.onequtritsquarepulse(ψI3, ν, Ω₀, T, ω₀, δ)

# FOR CHECKING ACCURACY
fidelity(ψ,φ) = 1 - abs2(ψ'*φ)

##########################################################################################
#                       PERFORM THE ACTUAL EVOLUTIONS - QUBIT CASE

println("Qubit Validation:")

# DIRECT EXPONENTIATION (Checked against ctrlq, this is the "standard".)
ψD = Evolutions.evolve(ψI2, pulses, device, Evolutions.Direct; numsteps=numsteps)
println("1-|⟨ψA2|ψD⟩|²: $(fidelity(ψA2, ψD))")

# ROTATE BETWEEN STATIC AND DRIVE BASES
ψRk = Evolutions.evolve(
    ψI2, pulses, device, Evolutions.Rotate; numsteps=numsteps,
    qubitapplymode=Evolutions.Kronec(), # APPLY QUBIT ROTATIONS WITH KRONECKER PRODUCT
)
println("1-|⟨ψA2|ψRk⟩|²: $(fidelity(ψA2, ψRk))")

ψRt = Evolutions.evolve(
    ψI2, pulses, device, Evolutions.Rotate; numsteps=numsteps,
    qubitapplymode=Evolutions.Tensor(), # APPLY QUBIT ROTATIONS WITH TENSOR ALGEBRA
)
println("1-|⟨ψA2|ψRt⟩|²: $(fidelity(ψA2, ψRt))")

# FACTOR DRIVE BASIS SO ALL ROTATIONS ARE TIME-INDEPENDENT
ψ1k = Evolutions.evolve(
    ψI2, pulses, device, Evolutions.Prediag; numsteps=numsteps,
    suzukiorder=1,                      # FACTOR WITH SIMPLEST POSSIBLE PRODUCT FORMULA
    qubitapplymode=Evolutions.Kronec(), # APPLY QUBIT ROTATIONS WITH KRONECKER PRODUCT
)
println("1-|⟨ψA2|ψ1k⟩|²: $(fidelity(ψA2, ψ1k))")

ψ1t = Evolutions.evolve(
    ψI2, pulses, device, Evolutions.Prediag; numsteps=numsteps,
    suzukiorder=1,                      # FACTOR WITH SIMPLEST POSSIBLE PRODUCT FORMULA
    qubitapplymode=Evolutions.Tensor(), # APPLY QUBIT ROTATIONS WITH TENSOR ALGEBRA
)
println("1-|⟨ψA2|ψ1t⟩|²: $(fidelity(ψA2, ψ1t))")

ψ2k = Evolutions.evolve(
    ψI2, pulses, device, Evolutions.Prediag; numsteps=numsteps,
    suzukiorder=2,                      # FACTOR WITH SYMMETRIC PRODUCT FORMULA
    qubitapplymode=Evolutions.Kronec(), # APPLY QUBIT ROTATIONS WITH KRONECKER PRODUCT
)
println("1-|⟨ψA2|ψ2k⟩|²: $(fidelity(ψA2, ψ2k))")

ψ2t = Evolutions.evolve(
    ψI2, pulses, device, Evolutions.Prediag; numsteps=numsteps,
    suzukiorder=2,                      # FACTOR WITH SYMMETRIC PRODUCT FORMULA
    qubitapplymode=Evolutions.Tensor(), # APPLY QUBIT ROTATIONS WITH TENSOR ALGEBRA
)
println("1-|⟨ψA2|ψ2t⟩|²: $(fidelity(ψA2, ψ2t))")

##########################################################################################
#                       PERFORM THE ACTUAL EVOLUTIONS - QUBIT CASE

println("Qutrit Validation:")

# DIRECT EXPONENTIATION (Checked against ctrlq, this is the "standard".)
ψD = Evolutions.evolve(ψI3, pulses, device, Evolutions.Direct; numsteps=numsteps)
println("1-|⟨ψA3|ψD⟩|²: $(fidelity(ψA3, ψD))")

# ROTATE BETWEEN STATIC AND DRIVE BASES
ψRk = Evolutions.evolve(
    ψI3, pulses, device, Evolutions.Rotate; numsteps=numsteps,
    qubitapplymode=Evolutions.Kronec(), # APPLY QUBIT ROTATIONS WITH KRONECKER PRODUCT
)
println("1-|⟨ψA3|ψRk⟩|²: $(fidelity(ψA3, ψRk))")

ψRt = Evolutions.evolve(
    ψI3, pulses, device, Evolutions.Rotate; numsteps=numsteps,
    qubitapplymode=Evolutions.Tensor(), # APPLY QUBIT ROTATIONS WITH TENSOR ALGEBRA
)
println("1-|⟨ψA3|ψRt⟩|²: $(fidelity(ψA3, ψRt))")

# FACTOR DRIVE BASIS SO ALL ROTATIONS ARE TIME-INDEPENDENT
ψ1k = Evolutions.evolve(
    ψI3, pulses, device, Evolutions.Prediag; numsteps=numsteps,
    suzukiorder=1,                      # FACTOR WITH SIMPLEST POSSIBLE PRODUCT FORMULA
    qubitapplymode=Evolutions.Kronec(), # APPLY QUBIT ROTATIONS WITH KRONECKER PRODUCT
)
println("1-|⟨ψA3|ψ1k⟩|²: $(fidelity(ψA3, ψ1k))")

ψ1t = Evolutions.evolve(
    ψI3, pulses, device, Evolutions.Prediag; numsteps=numsteps,
    suzukiorder=1,                      # FACTOR WITH SIMPLEST POSSIBLE PRODUCT FORMULA
    qubitapplymode=Evolutions.Tensor(), # APPLY QUBIT ROTATIONS WITH TENSOR ALGEBRA
)
println("1-|⟨ψA3|ψ1t⟩|²: $(fidelity(ψA3, ψ1t))")

ψ2k = Evolutions.evolve(
    ψI3, pulses, device, Evolutions.Prediag; numsteps=numsteps,
    suzukiorder=2,                      # FACTOR WITH SYMMETRIC PRODUCT FORMULA
    qubitapplymode=Evolutions.Kronec(), # APPLY QUBIT ROTATIONS WITH KRONECKER PRODUCT
)
println("1-|⟨ψA3|ψ2k⟩|²: $(fidelity(ψA3, ψ2k))")

ψ2t = Evolutions.evolve(
    ψI3, pulses, device, Evolutions.Prediag; numsteps=numsteps,
    suzukiorder=2,                      # FACTOR WITH SYMMETRIC PRODUCT FORMULA
    qubitapplymode=Evolutions.Tensor(), # APPLY QUBIT ROTATIONS WITH TENSOR ALGEBRA
)
println("1-|⟨ψA3|ψ2t⟩|²: $(fidelity(ψA3, ψ2t))")
