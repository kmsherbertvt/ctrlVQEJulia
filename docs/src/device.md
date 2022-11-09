# Devices
```@meta
CurrentModule = Devices
```

This module provides a framework for defining quantum computing devices.
The idea of a `Device` is to contain all information about the _static Hamiltonian_, such that
    knowing everything about the device tells you exactly how your qubits will evolve in time,
    _as long as you don't touch them_.

The "touching" generally happens through some form of electromagnetic radiation,
    and is treated elsewhere.

```@docs
Device
length(device::Device)
static_hamiltonian
Transmon
```
