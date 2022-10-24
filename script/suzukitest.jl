#= We want to see time and error as a function of Trotter steps,
    for each Evolution method: Trotter, Lanczos(p=0,1,2)
   We'll take Lanczos(p=0) as "exact" so don't plot error on that one... ^_^
=#

include("../src/utils.jl")
include("../src/device.jl")
include("../src/pulse.jl")
include("../src/evolve.jl")

using BenchmarkTools
import ..Devices
import ..Pulses
import ..Evolutions

fidelity(ψ,φ) = 1 - abs2(ψ' * φ)


##########################################################################################
#                               EXPERIMENTAL PARAMETERS
#=
        I'm just using the same ones I've been using,
        which are inspired by ctrlq examples.
        We'll vary only `numsteps`.
=#

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
#                                       RUN DATA
#=
        We'll run each method once to get the error, and then we'll run the fabulous
        @benchmark macro to collect as many time statistics as fit in a reasonable window.
=#
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 1.0



k = 1:14
t = Array{BenchmarkTools.Trial}(undef, (4,length(k)))
ε = Array{Float64}(undef, (4,length(k)))

for i in eachindex(k)
    numsteps = 2^k[i]
    println("Benchmarking k=$(k[i]), numsteps=$(numsteps)")

    # ERRORS
    ψ = Evolutions.evolve(ψI, pulses, device,
        Evolutions.Lanczos; suzukiorder=0,
        numsteps=numsteps,
    )
    ε[1,i] = 0
    ε[2,i] = fidelity(ψ, Evolutions.evolve(ψI, pulses, device,
        Evolutions.Lanczos; suzukiorder=1,
        numsteps=numsteps,
    ))
    ε[3,i] = fidelity(ψ, Evolutions.evolve(ψI, pulses, device,
        Evolutions.Lanczos; suzukiorder=2,
        numsteps=numsteps,
    ))
    ε[4,i] = fidelity(ψ, Evolutions.evolve(ψI, pulses, device,
        Evolutions.Trotter;
        numsteps=numsteps,
    ))

    # TIMES
    t[1,i] = @benchmark Evolutions.evolve(ψI, pulses, device,
        Evolutions.Lanczos; suzukiorder=0,
        numsteps=$numsteps,
    )
    t[2,i] = @benchmark Evolutions.evolve(ψI, pulses, device,
        Evolutions.Lanczos; suzukiorder=1,
        numsteps=$numsteps,
    )
    t[3,i] = @benchmark Evolutions.evolve(ψI, pulses, device,
        Evolutions.Lanczos; suzukiorder=2,
        numsteps=$numsteps,
    )
    t[4,i] = @benchmark Evolutions.evolve(ψI, pulses, device,
        Evolutions.Trotter;
        numsteps=$numsteps,
    )
end


#=
    Sometimes when ε gets very close to zero it can actually be negative.
    This indicates our states don't quite preserve normalization,
        which isn't great. That should maybe be addressed.
    But for now, just note that any ε less than the biggest occurrence
        should be taken with a grain of salt.
=#
ε_min = abs(minimum(ε))
ε .= abs.(ε)


# EXTRACT IMPORTANT STATISTICS FROM BENCHMARKING DATA
t_min = time.(minimum.(t))
t_med = time.( median.(t))


##########################################################################################
#                                       PLOT DATA
using Plots

labels = [
    "Lanczos p=0",
    "Lanczos p=1",
    "Lanczos p=2",
    "Trotter",
]
colors = [
    :purple,
    :blue,
    :green,
    :orange,
]

log10numsteps = log10.(2 .^ k)

# ERROR PLOT
errrplt = plot(title="Error", yaxis="log10| 1-|⟨ψL0|ψ⟩|² |", xaxis="log10 numsteps")
for i in 2:4
    plot!(log10numsteps, log10.(ε[i,:]), label=labels[i], color=colors[i])
end
hline!([log10(ε_min)], label=nothing)

# TIME PLOT
timeplt = plot(title="Time",yaxis="log10 t (ns)",xaxis="log10 numsteps",legend=:bottomright)
for i in 1:4
    plot!(log10numsteps, log10.(t_min[i,:]), label=labels[i], color=colors[i])
end
for i in 1:4
    plot!(log10numsteps, log10.(t_med[i,:]), label=nothing, color=colors[i], line=:dot)
end
