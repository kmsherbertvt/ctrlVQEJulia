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

Presently, this code features several methods to evolve a wave-function in time given a specific set of pulses,
but it does not yet implement any cost-function, gradient, or pulse optimization routines.
They are coming!
But in the meantime, consider Oinam Meitei's original [`ctrlq`](https://github.com/oimeitei/ctrlq) repository.

```@contents
Pages = ["utils.md", "device.md", "pulse.md", "evolve.md"]
```
