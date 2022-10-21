# ctrlVQEJulia

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





## Scripts
Scripts are located in the `script/` directory and are run from the project directory by ```julia script/<SCRIPT NAME>.jl>```.
- `trottertest` Consider a single time-step of the control Hamiltonian $H_C=\sum_q Ω_q(t) (e^{iν_qt} a_q + e^{-iν_qt} a_q^†)$.
  In principle, you should be able to approximate $\exp(-iΔt H_C)$ with a product formula, as long as $Δt$ is small enough.
  How small does it need to be? This script gives a brief sense of that, plotting a vague measure of error against $Δt$
  for two different product formulae. But really this was just my trial-by-fire for scientific computing in Julia...
- `evolvetest` I've replicated Oinam's evolution code; this script runs it and checks that its results more-or-less match.
  You can run the code yourself to see how much more-or-less is; it gets worse as you increase the duration `T`,
  but it's close enough that I think I'm okay attributing it to different round-off errors. Maybe.







## Usage Tutorial

### Install
```
> git clone https://github.com/kmsherbertvt/ctrlVQEJulia.git
> cd ctrlVQEJulia
```

I'm afraid that's it. I don't have the foggiest idea yet how to properly package Julia code, and I don't plan on having it anytime soon. ^_^

### Run Tests
```
> julia test/basicevolution.jl
```

This will run a few very basic unit tests, although it doesn't itself actually run any evolutions yet...

### Generate Documentation
If you have not installed Julia's `Documenter` package, type `julia` to enter the REPL, then type a `]` character to enter Pkg mode.
```
pkg> add Documenter
```
After it has installed, `<backspace>` gets you back to the regular REPL mode.

Now generate documentation via:
```
> cd docs
> julia make.jl
> cd ..
```

You can view documentation by opening the newly-created file `docs/build/index.html`.
Alas, browsers don't know to treat links to _local_ folders as opening their `index.html` files, so navigating the docs locally is a bit awkward.
The `Documenter` documentation offers some solutions to this (see the big blue note at the bottom of [this section](https://juliadocs.github.io/Documenter.jl/stable/man/guide/#Building-an-Empty-Document)), but personally I prefer to just manually open each index for now.

### Run Scripts
```
> julia script/<SCRIPT NAME>.jl>
```

In order to run `trottertest.jl`, you'll need the `Plots` package.
Type `julia` to enter the REPL, then type a `]` character to enter Pkg mode.
```
pkg> add Plots
```
After it has installed, `<backspace>` gets you back to the regular REPL mode.
