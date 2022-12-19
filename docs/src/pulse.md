# Pulses
```@meta
CurrentModule = Pulses
```

This module provides a framework for defining time-dependent pulses on a qubit.
The goal of the interface is to permit evolution and gradient algorithms
    to _not care_ how the time-dependent pulse is actually calculated.
You should not hesitate to define your own struct extending `PulseTemplate`
    to accommodate whatever parameterizations you are studying.

In the near future, the framework will be extended to include
    parametric pulse construction, amplitude, and frequency gradients.

```@docs
PulseTemplate
length(pt::PulseTemplate)
frequency
amplitude
BasicSquarePulse
```
