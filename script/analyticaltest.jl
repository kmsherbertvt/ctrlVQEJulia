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


println("\t1-|⟨ψA2|ψA2⟩|²: $(fidelity(ψA2, ψA2))")
println("\t1-|⟨ψA3|ψA3⟩|²: $(fidelity(ψA3, ψA3))")

##########################################################################################
#                       CONVENIENCE ANALYSIS FUNCTIONS

# FOR CHECKING ACCURACY
fidelity(ψ,φ) = 1 - abs2(ψ'*φ)

# FOR DISPLAYING RESULTS
function report(ψ, ψ0, header)
    println(header)
    println("\t1-|⟨ψ0|ψ⟩|²: $(fidelity(ψ0, ψ))")

    # # FOR CHECKING EFFECT OF NORMALIZATION
    # println()
    # println("\t1-|⟨ψ|ψ⟩|²: $(fidelity(ψ, ψ))")
    # println("\tForced normalization:")
    # ψ, ψ0 = normalize(ψ), normalize(ψ0)
    # println("\t1-|⟨ψ|ψ⟩|²: $(fidelity(ψ, ψ))")
    # println("\t1-|⟨ψ0|ψ⟩|²: $(fidelity(ψ0, ψ))")
    # println()
end

##########################################################################################
#                       PERFORM THE ACTUAL EVOLUTIONS - QUBIT CASE

println("Qubit Validation:")
println()

# DIRECT EXPONENTIATION (Checked against ctrlq, this is the "standard".)
ψD = Evolutions.evolve(ψI2, pulses, device, Evolutions.Direct; numsteps=numsteps)
report(ψD, ψA2, "DIRECT EXPONENTIATION")

# NUMERICAL INTEGRATION
ψS = Evolutions.evolve(ψI2, pulses, device, Evolutions.ODE;)
report(ψS, ψA2, "NUMERICAL INTEGRATION")

# ROTATE BETWEEN STATIC AND DRIVE BASES
ψRk = Evolutions.evolve(
    ψI2, pulses, device, Evolutions.Rotate; numsteps=numsteps,
    qubitapplymode=Evolutions.Kronec(), # APPLY QUBIT ROTATIONS WITH KRONECKER PRODUCT
)
report(ψRk, ψA2, "ROTATION (kronecker product)")

ψRt = Evolutions.evolve(
    ψI2, pulses, device, Evolutions.Rotate; numsteps=numsteps,
    qubitapplymode=Evolutions.Tensor(), # APPLY QUBIT ROTATIONS WITH TENSOR ALGEBRA
)
report(ψRt, ψA2, "ROTATION (tensor")

# FACTOR DRIVE BASIS SO ALL ROTATIONS ARE TIME-INDEPENDENT
ψ1k = Evolutions.evolve(
    ψI2, pulses, device, Evolutions.Prediag; numsteps=numsteps,
    suzukiorder=1,                      # FACTOR WITH SIMPLEST POSSIBLE PRODUCT FORMULA
    qubitapplymode=Evolutions.Kronec(), # APPLY QUBIT ROTATIONS WITH KRONECKER PRODUCT
)
report(ψ1k, ψA2, "PREDIAGONALIZED (order 1, kronecker product)")

ψ1t = Evolutions.evolve(
    ψI2, pulses, device, Evolutions.Prediag; numsteps=numsteps,
    suzukiorder=1,                      # FACTOR WITH SIMPLEST POSSIBLE PRODUCT FORMULA
    qubitapplymode=Evolutions.Tensor(), # APPLY QUBIT ROTATIONS WITH TENSOR ALGEBRA
)
report(ψ1t, ψA2, "PREDIAGONALIZED (order 1, tensor)")

ψ2k = Evolutions.evolve(
    ψI2, pulses, device, Evolutions.Prediag; numsteps=numsteps,
    suzukiorder=2,                      # FACTOR WITH SYMMETRIC PRODUCT FORMULA
    qubitapplymode=Evolutions.Kronec(), # APPLY QUBIT ROTATIONS WITH KRONECKER PRODUCT
)
report(ψ2k, ψA2, "PREDIAGONALIZED (order 2, kronecker product)")

ψ2t = Evolutions.evolve(
    ψI2, pulses, device, Evolutions.Prediag; numsteps=numsteps,
    suzukiorder=2,                      # FACTOR WITH SYMMETRIC PRODUCT FORMULA
    qubitapplymode=Evolutions.Tensor(), # APPLY QUBIT ROTATIONS WITH TENSOR ALGEBRA
)
report(ψ2t, ψA2, "PREDIAGONALIZED (order 2, tensor)")

##########################################################################################
#                       PERFORM THE ACTUAL EVOLUTIONS - QUBIT CASE

println("Qutrit Validation:")
println()

# DIRECT EXPONENTIATION (Checked against ctrlq, this is the "standard".)
ψD = Evolutions.evolve(ψI3, pulses, device, Evolutions.Direct; numsteps=numsteps)
report(ψD, ψA3, "DIRECT EXPONENTIATION")

# NUMERICAL INTEGRATION
ψS = Evolutions.evolve(ψI3, pulses, device, Evolutions.ODE;)
report(ψS, ψA3, "NUMERICAL INTEGRATION")

# ROTATE BETWEEN STATIC AND DRIVE BASES
ψRk = Evolutions.evolve(
    ψI3, pulses, device, Evolutions.Rotate; numsteps=numsteps,
    qubitapplymode=Evolutions.Kronec(), # APPLY QUBIT ROTATIONS WITH KRONECKER PRODUCT
)
report(ψRk, ψA3, "ROTATION (kronecker product)")

ψRt = Evolutions.evolve(
    ψI3, pulses, device, Evolutions.Rotate; numsteps=numsteps,
    qubitapplymode=Evolutions.Tensor(), # APPLY QUBIT ROTATIONS WITH TENSOR ALGEBRA
)
report(ψRt, ψA3, "ROTATION (tensor")

# FACTOR DRIVE BASIS SO ALL ROTATIONS ARE TIME-INDEPENDENT
ψ1k = Evolutions.evolve(
    ψI3, pulses, device, Evolutions.Prediag; numsteps=numsteps,
    suzukiorder=1,                      # FACTOR WITH SIMPLEST POSSIBLE PRODUCT FORMULA
    qubitapplymode=Evolutions.Kronec(), # APPLY QUBIT ROTATIONS WITH KRONECKER PRODUCT
)
report(ψ1k, ψA3, "PREDIAGONALIZED (order 1, kronecker product)")

ψ1t = Evolutions.evolve(
    ψI3, pulses, device, Evolutions.Prediag; numsteps=numsteps,
    suzukiorder=1,                      # FACTOR WITH SIMPLEST POSSIBLE PRODUCT FORMULA
    qubitapplymode=Evolutions.Tensor(), # APPLY QUBIT ROTATIONS WITH TENSOR ALGEBRA
)
report(ψ1t, ψA3, "PREDIAGONALIZED (order 1, tensor)")

ψ2k = Evolutions.evolve(
    ψI3, pulses, device, Evolutions.Prediag; numsteps=numsteps,
    suzukiorder=2,                      # FACTOR WITH SYMMETRIC PRODUCT FORMULA
    qubitapplymode=Evolutions.Kronec(), # APPLY QUBIT ROTATIONS WITH KRONECKER PRODUCT
)
report(ψ2k, ψA3, "PREDIAGONALIZED (order 2, kronecker product)")

ψ2t = Evolutions.evolve(
    ψI3, pulses, device, Evolutions.Prediag; numsteps=numsteps,
    suzukiorder=2,                      # FACTOR WITH SYMMETRIC PRODUCT FORMULA
    qubitapplymode=Evolutions.Tensor(), # APPLY QUBIT ROTATIONS WITH TENSOR ALGEBRA
)
report(ψ2t, ψA3, "PREDIAGONALIZED (order 2, tensor)")
