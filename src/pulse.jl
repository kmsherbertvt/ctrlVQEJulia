#= Implements different strategies to parameterize pulses Ω(t). =#

module Pulses

import ..Utils

"""
An abstract representation of the microwave pulse applied to a qubit.

A quantum control experiment will typically use an Array of PulseTemplates,
    with one PulseTemplate for each qubit.

# Implementation
A valid PulseTemplate `pt` should implement three methods:
- `Base.length(pt)` duration of pulse (alternatively, just include a `duration` attribute)
- `frequency(pt, t)` frequency of pulse at time `t`
- `amplitude(pt, t)` amplitude of pulse at time `t`

"""
abstract type PulseTemplate end


"""
    Base.length(pt::PulseTemplate)

Fetch the total time duration (ns) of this pulse.
"""
Base.length(pt::PulseTemplate) = pt.duration

"""
    frequency(pt::PulseTemplate, t::Float64)

Fetch the drive frequency (GHz) for this pulse at a specific time `t`.
"""
frequency(pt::PulseTemplate, t::Float64) = error("Not Implemented")

"""
    amplitude(pt::PulseTemplate, t::Float64)

Fetch the amplitude (GHz) for this pulse at a specific time `t`.

"""
amplitude(pt::PulseTemplate, t::Float64) = error("Not Implemented")

"""
    BasicSquarePulse(
        duration::Float64,
        frequency::Float64,
        amplitudes::Vector{Float64},
        steptimes::Vector{Float64},
    )

A basic square pulse has a constant frequency and amplitudes follow a step function.

# Arguments
- `duration` time duration (ns)
- `frequency` drive frequency (GHz)
- `amplitudes` list of `n` amplitudes (GHz)
- `steptimes` list of `n-1` times (ns) breaking up each amplitude, in ascending order

Times numerically equal to a breakpoint in `steptimes` are treated as after the breakpoint.

"""
struct BasicSquarePulse <: PulseTemplate
    duration::Float64
    frequency::Float64
    amplitudes::Vector{Float64}
    steptimes::Vector{Float64}

    # INNER CONSTRUCTOR: Validate lengths of `amplitudes` and `steptimes` vectors.
    BasicSquarePulse(duration, frequency, amplitudes, steptimes) = (
        length(steptimes) == length(amplitudes) - 1
            ? new(duration, frequency, amplitudes, steptimes)
            : error("BasicSqurePulse requires exactly one more amplitude than steptime.")
    )
end

"""
    BasicSquarePulse(duration::Float64, frequency::Float64, amplitude::Float64)

Construct a single-window square pulse.

"""
BasicSquarePulse(duration::Float64, frequency::Float64, amplitude::Float64) = (
    BasicSquarePulse(duration, frequency, [amplitude], Float64[])
)

frequency(pt::BasicSquarePulse, t::Float64) = pt.frequency

# amplitude(pt::BasicSquarePulse, t::Float64) = (
#     # DETERMINE WHICH AMPLITUDE TO TAKE, BASED ON WHERE `t` FALLS IN `timesteps`
#     ix = findfirst(steptime -> t < steptime, pt.steptimes);
#     # IF `ix` IS `nothing`, `t` IS AFTER ALL `timesteps`
#     return isnothing(ix) ? pt.amplitudes[end] : pt.amplitudes[ix]
# )

amplitude(pt::BasicSquarePulse, t::Float64) = sum(
    pt.amplitudes[i] * Utils.interval(t,
        get(pt.steptimes, i-1, 0),          #  LEFT TIME WINDOW (s[0] is implicitly 0)
        get(pt.steptimes, i, pt.duration),  # RIGHT TIME WINDOW (s[W] is implicitly T)
    ) for i in eachindex(pt.amplitudes)
)

end # END MODULE
