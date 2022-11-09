# Evolution
```@meta
CurrentModule = Evolutions
```

This module contains methods for evolving a quantum state in time,
    subject to a particular static Hamiltonian (encapsulated by a `device`)
    and a drive Hamiltonian (encapsulated by a list of `pulses`).

The mathematics gets dense here, so study the documentation and the code carefully,
    and feel free to ask for help if you feel like you don't understand something important!

```@docs
evolve!
evolve
```

Note that some additional auxiliary methods are also implemented in this module.
Study the documentation in code if you are curious,
    but they aren't necessarily meant for "public consumption". ^_^
