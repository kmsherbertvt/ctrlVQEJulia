#= We want to validate 'lanczos' evolution method gives same results as 'trotter'. =#

include("../src/utils.jl")
include("../src/device.jl")
include("../src/pulse.jl")
include("../src/evolve.jl")

module LanczosValidation
import ..Devices
import ..Pulses
import ..Evolutions

fidelity(ψ,φ) = 1 - abs2(ψ'*φ)

# CONSTRUCT INITIAL STATE
ψI::Vector{ComplexF64} = [0, 0, 0, 0, 1, 0, 0, 0, 0]

# DESIGN PULSE PATTERNS
T = 10.0
pulses = [
    Pulses.BasicSquarePulse(T, 29.0, [0.5, 0.4, 0.3, 0.2], [2.3, 5.6, 6.1]),
    Pulses.BasicSquarePulse(T, 31.0, [0.3, 0.1, 0.5, 0.6], [1.0, 3.4, 5.6]),
]

# DESIGN DEVICE
device = Devices.Transmon(
    2π .* [4.8080490154634950, 4.8332548172546130],
    2π .* [0.3101773613134229, 0.2916170385725456],
    Dict(
        Devices.QubitCouple(1,2) => 2π*0.018312874435769682,
    ),
)

##########################################################################################
#                           PERFORM THE ACTUAL EVOLUTIONS
ψT = Evolutions.evolve(ψI, pulses, device, Evolutions.Trotter; numsteps=2000)
ψL = Evolutions.evolve(ψI, pulses, device, Evolutions.Lanczos; numsteps=2000, suzukiorder=0)
ψ1 = Evolutions.evolve(ψI, pulses, device, Evolutions.Lanczos; numsteps=2000, suzukiorder=1)
ψ2 = Evolutions.evolve(ψI, pulses, device, Evolutions.Lanczos; numsteps=2000, suzukiorder=2)




##########################################################################################
#                           CHECK FIDELITY WITH EXPECTED RESULT

println("1-|⟨ψT|ψI⟩|²: $(fidelity(ψT, ψI))")
println("1-|⟨ψT|ψL⟩|²: $(fidelity(ψT, ψL))")
println("1-|⟨ψT|ψ1⟩|²: $(fidelity(ψT, ψ1))")
println("1-|⟨ψT|ψ2⟩|²: $(fidelity(ψT, ψ2))")

end
