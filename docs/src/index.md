# ctrl-VQE with Julia, m1

Welcome!
This code-base is meant to (one day) be a highly resource-efficient and flexible implementation of ctrl-VQE.

I put a very high value on thorough documentation,
via both this reference documentation and explanatory comments throughout the code.
If you find anything unclear (eg.
    "What on Earth is this code block doing?",
    "What does this variable name signify?",
    "Why do you calculate things this way?"
), please feel free to complain or file a Github issue:
I want this code to be accessible and understandable to anyone putting it to use.

The current code-base has a pretty well-tested evolution method,
    and a working amplitude gradient signal method.

The file `script/ctrlvqe_proofofconcept.jl` contains an in-depth tutorial on how to use these methods
    to run a VQE experiment with square pulses.

```@contents
Pages = ["utils.md", "device.md", "pulse.md", "evolve.md", "gradients.md"]
```
