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
```
