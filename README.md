# ctrlVQEJulia

Welcome!
This code-base is meant to (one day) be a highly resource-efficient and flexible implementation of ctrl-VQE.

I put a very high value on thorough documentation,
via both reference documentation (see below) and explanatory comments throughout the code.
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


## Dependencies

Naturally, the Julia language is required to run the code contained in the `ctrlVQEJulia` repository.

There are a whole bunch of packages used somewhere or another in the code.
I've lost track of them...
In the near future, this code base will be set up as a proper package,
  with all dependencies listed in a configuration file.

In the meantime, Julia will recognize all of the dependencies as available packages
  and offer to install them for you when you attempt to run the dependent code.

## Usage Tutorial

### Install
```
> git clone https://github.com/kmsherbertvt/ctrlVQEJulia.git
> cd ctrlVQEJulia
```

I'm afraid that's it.
I don't have the foggiest idea yet how to properly package Julia code,
  and I don't plan on having it anytime soon. ^_^

### Generate Documentation
HEY YOU YEAH YOU - I put a lot of effort into my documentation,
  so you should go through the steps here to build it,
  and then read it!
It should really help you "get" how this code works,
  and if it doesn't, I need you to tell me so I can rewrite it!

If you have not installed Julia's `Documenter` package, type `julia` to enter the REPL, then type a `]` character to enter Pkg mode.
```
pkg> add Documenter
```
After it has installed, `<backspace>` gets you back to the regular REPL mode. Type `exit()` to go back to shell.

Now generate documentation via:
```
> julia docs/make.jl
```

You can view documentation by opening the newly-created file `docs/build/index.html`.
Alas, browsers don't know to treat links to _local_ folders as opening their `index.html` files, so navigating the docs locally is a bit awkward.
The `Documenter` documentation offers some solutions to this (see the big blue note at the bottom of [this section](https://juliadocs.github.io/Documenter.jl/stable/man/guide/#Building-an-Empty-Document)), but personally I prefer to just manually open each index for now.

### Run Tests
```
> julia test/basicevolution.jl
```

This will run a few very basic unit tests, although it doesn't itself actually run any evolutions yet...

A proper software development framework writes tests _first_, and builds software around them.
I'm not doing that...

But, I do plan to have a good testing suite _some_ day.

### Run Scripts
```
> julia script/<SCRIPT NAME>.jl>
```

I create and discard scripts at will,
and generally make little effort to keep old scripts up-to-date with changes to code.
So, they are provided as-is with no expectation of functionality.

Nevertheless, they too are generally well-documented in the code,
so they may be interesting to glance at.
