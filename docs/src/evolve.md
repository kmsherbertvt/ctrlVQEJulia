# Evolution
```@meta
CurrentModule = Evolutions
```

```@docs
evolve(
    ψI::AbstractVector{<:Number},
    pulses::AbstractVector{<:Pulses.PulseTemplate},
    device::Devices.Device;
    numsteps::Integer = 500
)
evolve(
    ψI::AbstractVector{<:Number},
    pulses::AbstractVector{<:Pulses.PulseTemplate},
    device::Devices.Device,
    ::Type{Trotter};
    numsteps::Integer = 2000
)
evolve(
    ψI::AbstractVector{<:Number},
    pulses::AbstractVector{<:Pulses.PulseTemplate},
    device::Devices.Device,
    ::Type{Lanczos};
    numsteps::Integer = 2000,
    suzukiorder::Integer = 2,
)
```
