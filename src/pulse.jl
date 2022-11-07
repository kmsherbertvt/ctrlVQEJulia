#= Implements different strategies to parameterize pulses Ω(t). =#

module Pulses

"""
An abstract representation of the microwave pulse applied to a qubit.

A quantum control experiment will typically use an Array of PulseTemplates,
    with one PulseTemplate for each qubit.

# Implementation
A valid PulseTemplate pt should implement three methods:
- `Base.length(pt)`: duration of pulse (alternatively, just include a `duration` attribute)
- `frequency(pt, t)`: frequency of pulse at time `t`
- `amplitude(pt, t)`: amplitude of pulse at time `t`

"""
abstract type PulseTemplate end


"""
    Base.length(pt::PulseTemplate)

Fetch the total time duration (ns) of this pulse.
"""
Base.length(pt::PulseTemplate) = pt.duration

"""
    frequency(pt::PulseTemplate, t::Number)

Fetch the drive frequency (GHz) for this pulse at a specific time `t`.
"""
frequency(pt::PulseTemplate, t::Number) = error("Not Implemented")

"""
    amplitude(pt::PulseTemplate, t::Number)

Fetch the amplitude (GHz) for this pulse at a specific time `t`.

# TODO: How is an amplitude in GHz?
"""
amplitude(pt::PulseTemplate, t::Number) = error("Not Implemented")

"""
    BasicSquarePulse(
        duration::Number,
        frequency::Number,
        amplitudes::Vector{<:Number},
        steptimes::Vector{<:Number},
    )

A basic square pulse has a constant frequency and amplitudes follow a step function.

# Arguments
- `duration` time duration (ns)
- `frequency` drive frequency (GHz)
- `amplitudes` list of `n` amplitudes (GHz)
- `steptimes` list of `n-1` times (ns) breaking up each amplitude, in ascending order

Times numerically equal to a breakpoint in `steptimes` are treated as after the breakpoint.

# TODO: How is an amplitude in GHz? Or am I getting wires crossed...
"""
struct BasicSquarePulse <: PulseTemplate
    duration::Number
    frequency::Number
    amplitudes::AbstractVector{<:Number}
    steptimes::AbstractVector{<:Number}

    # INNER CONSTRUCTOR: Validate lengths of `amplitudes` and `steptimes` vectors.
    BasicSquarePulse(duration, frequency, amplitudes, steptimes) = (
        length(steptimes) == length(amplitudes) - 1
            ? new(duration, frequency, amplitudes, steptimes)
            : error("BasicSqurePulse requires exactly one more amplitude than steptime.")
    )
end

frequency(pt::BasicSquarePulse, t::Number) = pt.frequency

amplitude(pt::BasicSquarePulse, t::Number) = (
    # DETERMINE WHICH AMPLITUDE TO TAKE, BASED ON WHERE `t` FALLS IN `timesteps`
    ix = findfirst(steptime -> t < steptime, pt.steptimes);
    # IF `ix` IS `nothing`, `t` IS AFTER ALL `timesteps`
    return isnothing(ix) ? pt.amplitudes[end] : pt.amplitudes[ix]
)

end
