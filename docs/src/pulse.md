# Pulses
```@meta
CurrentModule = Pulses
```

This module provides a framework for defining time-dependent pulses on a qubit.
The goal of the interface is to permit evolution and gradient algorithms
    to _not care_ how the time-dependent pulse is actually calculated.
You should not hesitate to define your own struct extending `PulseTemplate`
    to accommodate whatever parameterizations you are studying.

TODO: haven't actually thought too much yet how to interface with a gradient algorithm.
It might be worth including "gradient" methods for amplitude and frequency some day.
They could default to a stock finite difference,
    but researchers can override the methods for their particular implementations
    to get more analytical results.

```@docs
PulseTemplate
length(pt::PulseTemplate)
frequency
amplitude
BasicSquarePulse
```
