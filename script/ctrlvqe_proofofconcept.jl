##########################################################################################
#= Welcome!

This file is meant to guide you through a sample ctrl-VQE calculation,
    using my Julia implementation.

The present version is somewhat transitory - it's not even a real package yet.
But since the code as of writing is *capable* of running ctrl-VQE
    (and much more efficiently than ctrlq) I wanted to get this version out
    to those of you doing *real* science (not just glorified software development ^_^).

=#########################################################################################


##########################################################################################
#                               PRELIMINARY SETUP
#=  The first thing to do is some generic Julia stuff, imports and such.                =#

# `using` STATEMENTS
using Random, LinearAlgebra         # STANDARD LIBRARIES
using Plots, LBFGSB                 # NEED TO BE INSTALLED
    #= `using` adds all of the [exported] names (ie. variables and methods)
                in the given package to the current namespace.

    I don't actually like using it because it makes tracing code significantly harder,
        so that my own code doesn't actually "export" anything at all.
    But we'll go ahead and use `using` for the standard libraries in this script,
        since that seems to be the Julia style. =#

# `include` STATEMENTS
include("../src/utils.jl")
include("../src/device.jl")
include("../src/pulse.jl")
include("../src/evolve.jl")
include("../src/gradient.jl")
include("../chemistry/molecule.jl")
    #= `include` is not quite the same as `import`.
            We'll actually have to import the contents of these files next.

    `include` is like copying the given source file directly into this script.
    It *shouldn't* be necessary,
        because I *should* be packaging my source code into packages.
    But I haven't gotten around to learning how to do that yet,
        so this is how it works for now. ^_^ =#

# `import` STATEMENTS
import ..Utils, ..Devices, ..Pulses, ..Evolutions, ..Gradients      # CORE MODULES
import ..Molecules                                                  # CHEMISTRY "PLUG-IN"
    #= `import` adds the given package name itself into the current namespace.
                We can access the package's members via dot notation (eg. `Utils.on`).

    The `..` preceding each package name is a little idiosyncracy due to the fact
        that these aren't actually packages (see `include` notes), but modules.
    That should be changed sometime, uh, soon-ish. =#

##########################################################################################
#                               CHEMISTRY SETUP
#=  Now let's define the problem. What system are we trying to solve?                   =#

#= TEMP  --- SKIP OVER THE pyscf PARTS

# CONSTRUCT A SECOND-QUANTIZED REPRESENTATION OF THE TARGET MOLECULE
geometry = Molecules.H2_geometry(0.7414)    # THE NUMBER IS THE ATOMIC SEPARATION
molecule = Molecules.Molecule(geometry)
    #= My `Molecules` module is a very slap-shod wrapper around pyscf,
        and I've only bothered to implement the bare minimum for H2.

    It is fairly straightforward to use for minimal-basis neutral singlet systems;
        you just need to specify the `geometry` corresponding to your system.
    You should be able to emulate the code in the `H2_geometry` function easily.

    For charged systems or those with uneven spins or larger basis sets, sorry -
        I was too lazy to include these parameters in my `Molecule` constructor. :P
    If you're familiar with pyscf, you can probably edit my code yourself to include them.
        There really shouldn't be anything to it, just adding a few lines of code.

    Otherwise, the easiest strategy is probably to skip the above step and
        load in your Hamiltonian matrix directly from an .npz file,
        produced from a Python script or obtained from the group.

    You can use Julia's `NPZ` package to load it in as a Julia matrix.
        But you'll have to look up the syntax yourself. ;) =#

# CONSTRUCT THE FERMIONIC HAMILTONIAN IN THE COMPUTATIONAL BASIS
H = Molecules.molecular_hamiltonian(molecule)
n = molecule.n

=#

# LOAD THE HAMILTONIAN FROM AN `.npy` FILE
using NPZ
H = npzread("matrix/h215.npy")
n = round(Int, log2(size(H,1)))


# INITIAL KET |1100⟩
init_ket = "1"^(n÷2) * "0"^(n÷2)
    #= Later on, we'll initialize a statevector |ψ0⟩ so that
        ⟨z|ψ0⟩=1 if |z⟩ is |init_ket⟩, and ⟨z|ψ0⟩=0 otherwise.

    Each character represents the occupation of a spin-orbital.
    Odd-indexed characters correspond to α spin; even-indexed to β spin.
        (NOTE: Does that sound backwards? Julia indexes from 1, not 0!)

    The usual thing is to set |init_ket⟩ to the Hartree Fock state,
        which fills up orbitals from the left with available electrons.
    For example, in a singlet system, there are as many α-electrons as β-electrons,
        so the Hartree-Fock ket is just as many 1's as there are electrons,
        followed by as many 0's as needed to address the remaining spin-orbitals. =#

##########################################################################################
#                                  PROBLEM SPECS
#=  Here we'll define a bunch of variables that will be used throughout.                =#

# THE NUMBER OF MODES CONSIDERED FOR EACH TRANSMON
m = 2
    #= Note the distinction:
        In a digital quantum computation, every qubit has two levels,
            and it's easy(ish) to map each qubit to the fermionic electron problem.
        But in ctrl-VQE, we are closer to hardware:
            we are forced to treat qubits as *bosonic* systems,
            with an infinite-dimensional Hilbert space.
        We have to truncate at some point: this is `m`.

    Note if our pulses never ever produce any leakage
            into modes outside the computational space (ie. 0 and 1 on each transmon),
            m=2 is sufficient, and easiest to simulate.
    BUT there could be some benefit to intentionally inducing leakage;
        that's one of our lines of research! =#

# PHYSICAL CHARACTERISTICS OF OUR DEVICE
device = Devices.selectqubits(1:n, Devices.Transmon(
    2π*[3.7, 4.2, 3.5, 4.0],                    # QUBIT RESONANCE FREQUENCIES
    2π*[0.3, 0.3, 0.3, 0.3],                    # QUBIT ANHARMONICITIES
    Dict{Devices.QubitCouple,Float64}(          # QUBIT COUPLING CONSTANTS
        Devices.QubitCouple(1,2) => 2π*.018,
        Devices.QubitCouple(2,3) => 2π*.021,
        Devices.QubitCouple(3,4) => 2π*.020,
        Devices.QubitCouple(1,3) => 2π*.021,
        Devices.QubitCouple(2,4) => 2π*.020,
        Devices.QubitCouple(1,4) => 2π*.021,
    )
))
    #= There are a whole bunch of things to note here:

    1. The numbers in the example are utterly arbitrary.
       Orders of magnitude are (probably) consistent with IBMQ devices,
        but the numbers were haphazardly generated at whim.
       You should absolutely change them to something sensible. ^_^

    1b. I inherited the convention of using units of 2π from ctrlq;
        I don't really understand why I'm doing it... ^_^

    2. The example is a device with four qubits.
       If your calculation needs more qubits,
        just add on to the frequency and anharmonicity vectors accordingly,
        and add qubit couplings as you deem appropriate.
       If your calculation needs fewer qubits,
        the call to `selectqubits` takes care of it.

    3. The qubit coupling constants are assumed to be symmetric.
       Thus, you could specify the pairing between qubits 1 and 2 with EITHER of
        `Devices.QubitCouple(1,2)`  OR  `Devices.QubitCouple(2,1)`
       They are the same object.

    4. Julia indexes from 1, not from 0.
        For example, `Devices.QubitCouple(1,2)` indicates the pairing between
            the *first and second* qubits. =#

# NUMBER OF TIME STEPS
r = 200
    #= In the ctrlq code, a typical number of time steps was on the order of 1000.
        This code has made some numerical improvements,
            so the same results should be obtainable with a fraction of that.

    But if something isn't working the way you expect,
        it's not a bad idea to bump this up a couple orders of magnitude
        and see if that fixes it... =#


# PULSE CONSTRAINTS
ΔΩ = 2π * 0.2      # RANGE OF PULSE AMPLITUDE
ΔΔ = 2π * 1.0       # RANGE OF DE-TUNING
    #= Pulse amplitude Ω on each qubit will be constrained to `Ω ∈ (-ΔΩ, +ΔΩ)`.

    Pulse frequency ν on qubit with resonance ω will be constrained to `ν ∈ (ω-ΔΔ, ω+ΔΔ)`.
        I don't know if constraining frequency like this makes sense in experiment,
            but since the pulse gets less and less *relevant* as Δ≡ω-ν gets larger,
            it makes sense in theory. =#


# PULSE HYPER-PARAMETERS
T = 20.0        # TOTAL PULSE DURATION
W = 40         # NUMBER OF PULSE WINDOWS
    #= You'll likely want to change these,
        or write your own script to change them adaptively.
       This is just a proof-of-concept, after all! =#
steptimes = repeat( range(0,T,W+1)[2:end-1], 1, n)      # PULSE WINDOW STEP TIMES
    #= This is an n-column matrix setting each window to be equally spaced.
        You can manually set it, or write your own adaptive strategy. =#



##########################################################################################
#                                  MORE SETUP!
#=  Now we're going to derive a bunch of things, from the parameters given above.
    You shouldn't need to change anything here,
        but study the notes so you know what's going on!                                =#

# TOTAL SIZE OF HIBLERT SPACE BEING SIMULATED
N = m^n
    #= NOTE this is generally larger than the 2^n Hilbert space being "solved". =#

# SINGLE-QUBIT BOSONIC ANNIHILATION OPERATOR (truncated to m modes)
a = Utils.a_matrix(m)

# EXTEND OBSERVABLE ONTO SIMULATION SPACE
Π = Utils.projector(n, 2, m)    # PROJECTOR FROM SIMULATION SPACE ONTO PROBLEM SPACE
O = Hermitian(Π'*H*Π)           # MOLECULAR HAMILTONIAN, EMBEDDED INTO SIMULATION SPACE
    #= Nick tells me that Π is really an "isometry", and the projector is Π'Π.
        Names will likely change in a future version.

    The `Hermitian` command just wraps the matrix Π'HΠ into a "Hermitian" matrix object,
        so that Julia knows it can use more efficient algorithms on it. =#

# PREPARE THE HARTREE-FOCK STATE
ψ0 = zeros(ComplexF64, N)                               # INITIALIZE ⟨z|ψ0⟩ = 0
ψ0[1 + parse(Int, init_ket, base=m)] = one(ComplexF64)  # SET ⟨init_ket|ψ0⟩ = 1
    #= Here we are using the `init_ket` set in the chemistry section,
        to initialize the corresponding statevector in the simulation space.

    The `parse(Int, init_ket, base=m)` interprets the `init_ket` string
        as an integer represented with base-`m` digits.
    We have to add `1 +` because Julia starts indexing from 1,
        so for example |000⟩ is actually index 1. =#

# SET OPTIONS FOR DIFFERENT CALCULATION MODES
qubitapplymode = Evolutions.Kronec()    # HOW TO APPLY SINGLE-QUBIT OPERATORS?
iobasis = Evolutions.QubitBasis()       # HOW TO INTERPRET THE INPUT/OUTPUT STATEVECTOR?
    #=
    `qubitapplymode` could be `Evolutions.Kronec()` or `Evolutions.Tensor()`

        In principle, `Tensor` mode is faster,
            but there's such significant overhead that up to ~5 qubits
            you probably want to stick with `Kronec`.

        In truth, single-qubit operations aren't the "rate-limiting step", so to speak,
            so when in doubt, stick with `Kronec`. ;)

    `iobasis` could be `Evolutions.QubitBasis()` or `Evolutions.DeviceBasis()`

        We're not entirely sure yet which one more closely matches a real experiment.

        It (no longer) makes any computational difference,
            so I've picked `QubitBasis` because I find it more intuitive. =#

# DISCRETIZE TIME
t_= range(0,T,r+1)                 # TIME GRID
τ = T / r                          # DURATION OF EACH TIME STEP
    #= I like to use an underscore suffix to represent an array.
        So, `t` is a single time point, and `t_` is an array of time points.

    Note that `τ == t_[2]-t[1]`, and there are a total of `r` "steps" in `t_`. =#

# HEAVY CALCULATIONS FOR WORKING IN THE DEVICE BASIS (aka "dressed basis")
HD = Devices.static_hamiltonian(device, m)  # DEVICE HAMILTONIAN
ΛD, UD = Utils.dressedbasis(HD)             # DEVICE EIGENVALUES AND EIGENVECTORS
L = UD* Diagonal(exp.((-im*τ) * ΛD)) *UD'   # REPEATED DEVICE ACTION FOR EACH TIME STEP
    #= The method `Utils.dressedbasis` is essentially just factoring H→UΛU',
        but there are a couple extra steps, like changing the ordering of eigenvectors.

    The letter `L` stands for "ligand", because `L` is the "ligand operator"
        connecting the application of the drive Hamiltonian at each time step. =#

# PRE-ALLOCATIONS
ψ  = Vector{ComplexF64}(undef, m^n)         # FOR STORING THE WAVEFUNCTION
∇Ω = Matrix{Float64}(undef, r+1, n)         # FOR STORING THE DERIVATIVES
    #= The important calculations will be written to these variables. =#

_N1 = Vector{ComplexF64}(undef, N)          # FOR MATRIX-VECTOR MULTIPLICATION
_m2_ = [Matrix{ComplexF64}(undef, m,m) for q ∈ 1:n]     # QUBIT-WISE DRIVE OPERATORS
_K_ = [Matrix{ComplexF64}(undef, m^q, m^q) for q ∈ 1:n] # FOR KRON'ing OPERATORS
    #= To minimize the computational effort at each time-step,
        the evolve and gradient functions let you pass "pre-allocated" blocks of memory
        that they can use as "work" variables, to store intermediate results in.

    I'm prefixing these with an underscore, to flag them as "meaningless to inspect".
    The names represent their use:
        `_N1` is a 1d vector of `N` elements.
        `_m2_` is an *array* of `m×m` matrices.
        `_K_` is the most niche; it's an array of matrices of increasing size.
            Each one is used to store an intermediate result in the
                (K)ronecker products of many single-qubit operators. =#

# WRAP ALL THESE UP INTO A NICE, TIDY KEYWORD DICTIONARY
evolve_kwargs = Dict(
    :iobasis => iobasis, :qubitapplymode => qubitapplymode, :numsteps => r,
    :N => N, :n => n, :m => m, :T => T, :t_=> t_, :Δt=> τ,
    :ΛD => ΛD, :UD => UD, :V  => L, :a  => a,
    :tmpV  => _N1, :tmpM_ => _m2_, :tmpK_ => _K_,
)

gradient_kwargs = Dict(
    :iobasis => iobasis, :qubitapplymode => qubitapplymode, :r => r,
    :N => N, :n => n, :m => m, :T => T, :t_=> t_, :τ => τ,
    :ΛD => ΛD, :UD => UD, :V  => L, :a  => a,
    :tmpV  => _N1, :tmpM_ => _m2_, :tmpK_ => _K_, :∂Ω => ∇Ω,
)
    #= The evolve and gradient functions will re-compute a lot of this stuff automatically
        unless you provide pre-calculated values as keywork arguments.
    But, providing all these keyword arguments makes the function calls look very ugly.
    So I'm wrapping it all into a nice dictionary
        that can just be dereferenced in the function calls. =#


##########################################################################################
#                    DEFINE PULSE, COST FUNCTION, AND GRADIENT FUNCTION
#=  My apologies - there is more happening in this section than I had intended.
    I had meant for many of these things to be handled inside my core modules,
        but I realized that the way I want to handle them will require big changes,
    and I want to get this tutorial out ASAP.

    Therefore I'm implementing the details of working with a square pulse here,
        and including lots of notes so that you understand
        how to adapt it to any pulse shape you like. ^_^

    In a later version, these things should all be bound up into a tidy module,
        which you can simply extend as you define new pulse shapes.                     =#


"""
    constructpulse(x::AbstractVector{<:Real})

This function constructs a set of pulse objects for each qubit,
    out of a vector of optimizable parameters.

In this example, each pulse will be a sequence of `W` square pulses.

The vector `x` must specify:
- `W` different amplitudes A[j]
- a single ν, for the frequency of the pulse.

Actually, `x` must specify each of these for all `n` qubits,
    for a total of `n(W+1)` parameters.

"""
function constructpulses(x::AbstractVector{<:Real})
    x = reshape(x, :, n)    # RE-INTERPRET `x` AS `n` COLUMNS

    pulses = Vector{Pulses.BasicSquarePulse}(undef, n)
    for q in 1:n
        A = x[1:W,q]                # AMPLITUDE ARRAY
        ν = x[end,q]                # FREQUENCY
        pulses[q] = Pulses.BasicSquarePulse(T, ν, A, steptimes[:,q])
    end
    return pulses
end


"""
    pulsegradient_amplitude(t, j, pulse)

Consider the pulse amplitude Ω(t) on a single qubit.
    In our example, Ω(t) is a step function, determined by each A[j] and s[j].
    We're assuming s[j] are held constant in this tutorial,
        but a bit of calculus gives us the partial derivatives ∂Ω(t)/∂A[j].

This function gives ∂Ω(t)/∂A[j].

"""
function pulsegradient_amplitude(t, j, pulse)
    return Utils.interval(t,
        get(pulse.steptimes, j-1, 0),               #  LEFT BOUND (TREAT s[0] as t=0)
        get(pulse.steptimes, j, pulse.duration),    # RIGHT BOUND (TREAT s[W] as t=T)
    )
end

∂Ω∂x = Vector{Float64}(undef, r+1)  # PRE-ALLOCATE A VECTOR TO HOLD THESE GRADIENTS
h = 1e-6                # SPACING FOR A FINITE DIFFERENCE ON THE FREQUENCY PARAMETER


"""
    f(x)

The cost function to be used in optimization.

Returns the energy ⟨ψ|O|ψ⟩, where |ψ⟩ is the quantum state prepared
    by applying the pulse sequences parameterized by `x` to the initial state |ψ0⟩.

Note that this energy is "unnormalized",
    meaning it will actually be implicitly a little higher
    than the one you'd get by ⟨ψ'|H|ψ'⟩, where |ψ'⟩ is |ψ⟩ retricted to 0's and 1's.
That means leakage into modes larger than 0 and 1 will be implicitly penalized
    in an optimization.

The true minimum of O will be exactly equal to the true minimum of H,
    but there's no guarantee an optimization *finds* the true minimum,
    so you must correct for leakage before interpreting the accuracy of your optimization:

    ⟨ψ'|H|ψ'⟩ = ⟨ψ|O|ψ⟩ / ⟨ψ|Π'Π|ψ⟩

"""
function f(x)
    pulses = constructpulses(x)
    ψ .= Evolutions.evolve(ψ0, pulses, device, Evolutions.Rotate; evolve_kwargs...)
    return real(Utils.expectation(O, ψ))
end


"""
    g!(G, x)

The gradient function to be used in optimization.

`G` is just a pre-allocated block of memory to write the gradient vector to.
The `!` in the name is just Julia convention to mark this function as one
    which mutates its primary operand (in this case, `G`).

To calculate the gradients ∂E/∂A[j],
    we'll need to do some rather unrigorous calculus:

    ∂E/∂x[j] = ∂E/∂Ω(t) · ∂Ω(t)/∂x[j]
                        ↑
            This is a vector dot product.

We've implemented ∂Ω(t)/∂x[j] at the start of this section.
The ∂E/∂Ω(t) is what we're calling the gradient signal, or switching function φ(t).

    (Actually my gradient signal has an extra `τ` (aka `dt`) in it,
        and I'm not exactly certain what it's doing there.
    But it seems to work... )

To calculate the gradients ∂E/∂ν...
    ah, well, there is a *frequency* gradient signal I could calculate.
But I haven't implemented that yet.
So what we're going to do instead is run a forward finite difference,
    on just these `n` parameters... ^_^

"""
function g!(G, x)
    G = reshape(G, :, n)    # RE-INTERPRET `G` AS `n` COLUMNS

    # CONSTRUCT THE AMPLITUDE GRADIENT SIGNAL ∇Ω ~= ∂E/∂Ω(t) dt
    pulses = constructpulses(x)
    Gradients.gradientsignal(ψ0, pulses, device, O; gradient_kwargs...)
        #= Because we put a pre-allocated `∇Ω` into `gradient_kwargs`,
            this method writes the gradient signal to that `∇Ω`. =#

    # EVALUATE THE COST FUNCTION AT x, FOR FORWARD FINITE DIFFERENCES
    f0 = f(x)

    for (q, pulse) in enumerate(pulses)
        for j in 1:W
            # SET THE GRADIENTS G[j] ~= ∂E/∂A[j]
            ∂Ω∂x .= map(t -> pulsegradient_amplitude(t,j,pulse), t_)
            G[j,q] = ∇Ω[:,q] · ∂Ω∂x
        end

        # SET THE GRADIENTS G[j] ~= ∂E/∂ν
        Δν = h * Utils.basisvector(length(x), (W+1)*q)
        G[W+1,q] = (f(x+Δν) - f0) / h
    end

end



##########################################################################################
#                               RUN L-BFGS-B OPTIMIZATION
#=  In this section, we'll set up our constraints and initial guess,
        and run the optimization.

    The optimization package I'm using is a wrapper around C code, I think. =#


# CONFIGURE BOUNDS
lbounds = Matrix{Float64}(undef, W+1, n); ubounds = Matrix{Float64}(undef, W+1, n)
# # #
lbounds[1:W,:] .= -ΔΩ                       # SET MAXIMUM AMPLITUDE
ubounds[1:W,:] .=  ΔΩ
# # #
lbounds[end,:] .= device.ω .- ΔΔ            # SET FREQUENCY CONSTRAINTS
ubounds[end,:] .= device.ω .+ ΔΔ


# CONFIGURE INITIAL PARAMETERS
Random.seed!(0)                         # CHANGE SEED AT WILL!
x0 = Matrix{Float64}(undef, W+1, n)
# # #
x0[1:W,:] .= ΔΩ * (2*rand(W,n) .- 1)                        # RANDOMLY SELECTED AMPLITDUES
x0[end,:] .= ΔΔ * (2*rand(n) .- 1) .+ device.ω              # RANDOMLY SELECTED FREQUENCIES
x0 = vec(x0)                            # RESHAPE TO VECTOR

# CALCULATE INITIAL ENERGY
E0 = f(x0)
println("     Initial energy: $E0")

E_HF = f(zero(x0))
println("Hartree Fock energy: $E_HF")

# PLOT PULSE SHAPES AND GRADIENT SIGNALS
pulses = constructpulses(x0)            # INITIAL PULSE OBJECTS
Ω0 = Matrix{Float64}(undef, r+1, n)     # CONSTRUCT Ω(t) FOR EACH QUBIT
for (q, pulse) in enumerate(pulses)
    Ω0[:,q] .= map(t -> Pulses.amplitude(pulse, t), t_)
end
# # #
Ω_plots = plot(                         # PULSE SHAPE PLOT
    [plot(
        t_, Ω0[:,q]
    ) for q in 1:n]...,
    title = "Initial Pulse Shapes",
    ylim = [-ΔΩ, ΔΩ],
    legend = false,
    layout = (n,1),
)
# # #
Gradients.gradientsignal(ψ0, pulses, device, O; gradient_kwargs...) # GRADIENT SIGNAL
∇Ω_plots = plot(                        # GRADIENT SIGNAL PLOT
    [plot(
        t_, ∇Ω[:,q]
    ) for q in 1:n]...,
    title = "Initial Gradient Signals",
    legend = false,
    layout = (n,1),
)
# # #
plot(Ω_plots, ∇Ω_plots, layout=(1,2))   # COMBINE PLOTS INTO ONE
gui()                                   # *DISPLAY* PLOTS

# RUN THE OPTIMIZATION
_, x = lbfgsb(f, g!, vec(x0), lb=vec(lbounds), ub=vec(ubounds),
# _, x = lbfgsb(f, ng!, vec(x0), lb=vec(lbounds), ub=vec(ubounds),
    m=5,            # use last m iterations to inform new α, or something...
    factr=1e1,      # ENERGY CONVERGENCE CRITERION, in units of machine epsilon
    pgtol=1e-5,     # GRADIENT CONVERGENCE CRITERION (norm of projected gradient)
    iprint=50,      # PRINT EVERY FEW ITERATIONS, but it means something different for 99+
    maxfun=10000,   # How many times to call f before we give up?
    maxiter=10000,  # How many linesearches to run before we give up?
)
    #= `lbfgsb`s first return value (suppressed with `_`) is the energy `Ex`,
            but by calling `f(x)` manually at the end,
            I guarantee that the last values written to `ψ` are the final statevector. =#


##########################################################################################
#                                    FINAL RESULTS



# CONSTRUCT A BUNCH OF USEFUL VALUES FROM OPTIMIZATION RESULTS
pulses = constructpulses(x)         # CONSTRUCT FINAL PULSES
Gradients.gradientsignal(ψ0, pulses, device, O; gradient_kwargs...) # + GRADIENT SIGNAL
Ex = f(x)                           # CONSTRUCT FINAL ENERGY AND STATEVECTOR

leakage = 1 - real(ψ'*Π'*Π*ψ)       # CALCULATE LEAKAGE
E_normalized = Ex / (1 - leakage)   # CALCULATE RE-NORMALIZED ENERGY

Ω = Matrix{Float64}(undef, r+1, n)  # CONSTRUCT Ω(t) FOR EACH QUBIT
for (q, pulse) in enumerate(pulses)
    Ω[:,q] .= map(t -> Pulses.amplitude(pulse, t), t_)
end

# DISPLAY OPTIMAL FREQUENCIES
println("="^30)
println("Optimized frequencies for each qubit:")
for (q, pulse) in enumerate(pulses)
    println("\tQubit $q: ν=$(pulse.frequency)")
end

# REPORT ENERGIES
println("="^30)
println("            Leakage: $leakage")
println("Unnormalized energy: $Ex")
println("Renormalized energy: $E_normalized")

# PLOT PULSE SHAPES AND GRADIENT SIGNALS
Ω_plots = plot(
    [plot(
        t_, Ω[:,q]
    ) for q in 1:n]...,
    title = "Pulse Shapes",
    ylim = [-ΔΩ, ΔΩ],
    legend = false,
    layout = (n,1),
)

∇Ω_plots = plot(
    [plot(
        t_, ∇Ω[:,q]
    ) for q in 1:n]...,
    title = "Gradient Signals",
    legend = false,
    layout = (n,1),
)

plot(Ω_plots, ∇Ω_plots, layout=(1,2))





# CONSTRUCT A COUPLE USEFUL VALUES FROM SOLVING FCI MATRIX (if tractable)
Λ, U = eigen(O)
E_FCI = Λ[1]
ψ_FCI = U[:,1]

println("                FCI energy: $E_FCI")
println("  Unoptimized energy error:  $(E0-E_FCI)")
println(" Hartree Fock energy error:  $(E_HF-E_FCI)")
println()
println("    Optimized energy error:  $(E_normalized-E_FCI)")
println("Optimized state infidelity:  $(Utils.infidelity(ψ,ψ_FCI))")
println()
println("       Energy error wrt HF:  $((E_normalized-E_FCI)/(E_HF-E_FCI))")