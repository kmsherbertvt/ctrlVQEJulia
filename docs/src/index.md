# ctrl-VQE with Julia, m1

Welcome! This is Kyle's first attempt at implementing ctrl-VQE in Julia.
The goal is to produce a working code giving results numerically similar to Oinam's `ctrlq` repository,
    and to test out a presumably-faster Lanczos method rather than classic Trotterization.

(They're actually the exact same method except the classic "Trotter" version operates in the device basis
    and treats each time step as independent, leading to rather more matrix exponentiations than are needed.
The "Lanczos" method will work in the computational basis, enabling somewhat easier treatment of the control Hamiltonian.
Meanwhile, the action of the device Hamiltonian is actually identical at each time step (so long as they're all the same length),
    so we only actually need to compute it once.)

As this is Kyle's first foray into Julia, there is no expectation of perfectly-optimized code.
In particular, exactly zero effort is being made to make this code parallelizable.
That will come in a lower-priority second attempt.

```@contents
```
