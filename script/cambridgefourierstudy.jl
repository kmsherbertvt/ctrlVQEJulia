#= Run optimization up to a point, and publish figures of time and frequency domain. =#

maxiter = 256     # WHEN TO TERMINATE OPTIMIZATION



T = 30.0        # TOTAL PULSE DURATION
W = 50         # NUMBER OF PULSE WINDOWS
r = 10000       # NUMBER OF TROTTER STEPS



systemtag = "lih30"
ΔΩ = 2π * .02       # RANGE OF PULSE AMPLITUDE
ΔΔ = 2π * 1.0       # RANGE OF DE-TUNING
m = 2           # NUMBER OF LEVELS PER TRANSMON
seed = 5        # RANDOM SEED TO SELECT INITIAL PARAMETERS

##########################################################################################

# `using` STATEMENTS
using Random, LinearAlgebra         # STANDARD LIBRARIES
using Plots, LBFGSB                 # NEED TO BE INSTALLED

# `include` STATEMENTS
include("../src/utils.jl")
include("../src/device.jl")
include("../src/pulse.jl")
include("../src/evolve.jl")
include("../src/gradient.jl")
include("../chemistry/molecule.jl")

# `import` STATEMENTS
import ..Utils, ..Devices, ..Pulses, ..Evolutions, ..Gradients      # CORE MODULES

##########################################################################################

# LOAD THE HAMILTONIAN FROM AN `.npy` FILE
using NPZ
H = npzread("matrix/$systemtag.npy")
E_FCI = eigen(H).values[1]
n = round(Int, log2(size(H,1)))

# PHYSICAL CHARACTERISTICS OF OUR DEVICE
device = n <= 4 ?
    Devices.selectqubits(1:n, Devices.NPJ_DEVICE) :
    Devices.SYSTEMATIC_DEVICE(n)

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

# PREPARE THE HARTREE-FOCK STATE
ψ0 = zeros(ComplexF64, N)                               # INITIALIZE ⟨z|ψ0⟩ = 0
ψ0[argmin(real.(diag(O)))] = one(ComplexF64)            #    ASSIGN ⟨HF|ψ0⟩ = 1

# SET OPTIONS FOR DIFFERENT CALCULATION MODES
qubitapplymode = Evolutions.Kronec()    # HOW TO APPLY SINGLE-QUBIT OPERATORS?
iobasis = Evolutions.QubitBasis()       # HOW TO INTERPRET THE INPUT/OUTPUT STATEVECTOR?

# DISCRETIZE TIME
t̄ = range(0,T,r+1)                 # TIME GRID
τ = T / r                          # DURATION OF EACH TIME STEP

# HEAVY CALCULATIONS FOR WORKING IN THE DEVICE BASIS (aka "dressed basis")
HD = Devices.static_hamiltonian(device, m)  # DEVICE HAMILTONIAN
ΛD, UD = Utils.dressedbasis(HD)             # DEVICE EIGENVALUES AND EIGENVECTORS
L = UD* Diagonal(exp.((-im*τ) * ΛD)) *UD'   # REPEATED DEVICE ACTION FOR EACH TIME STEP

# PRE-ALLOCATIONS
ψ  = Vector{ComplexF64}(undef, m^n)         # FOR STORING THE WAVEFUNCTION
∇Ω = Matrix{Float64}(undef, r+1, n)         # FOR STORING THE DERIVATIVES

_N1 = Vector{ComplexF64}(undef, N)          # FOR MATRIX-VECTOR MULTIPLICATION
_m2_ = [Matrix{ComplexF64}(undef, m,m) for q ∈ 1:n]     # QUBIT-WISE DRIVE OPERATORS
_K_ = [Matrix{ComplexF64}(undef, m^q, m^q) for q ∈ 1:n] # FOR KRON'ing OPERATORS

# WRAP ALL THESE UP INTO A NICE, TIDY KEYWORD DICTIONARY
evolve_kwargs = Dict(
    :iobasis => iobasis, :qubitapplymode => qubitapplymode, :numsteps => r,
    :N => N, :n => n, :m => m, :T => T, :t_ => t̄, :Δt=> τ,
    :ΛD => ΛD, :UD => UD, :V  => L, :a  => a,
    :tmpV  => _N1, :tmpM_ => _m2_, :tmpK_ => _K_,
)

gradient_kwargs = Dict(
    :iobasis => iobasis, :qubitapplymode => qubitapplymode, :r => r,
    :N => N, :n => n, :m => m, :T => T, :t_ => t̄, :τ => τ,
    :ΛD => ΛD, :UD => UD, :V  => L, :a  => a,
    :tmpV  => _N1, :tmpM_ => _m2_, :tmpK_ => _K_, :∂Ω => ∇Ω,
)

# PULSE WINDOW STEP TIMES
steptimes = repeat( range(0,T,W+1)[2:end-1], 1, n)


##########################################################################################


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
            ∂Ω∂x .= map(t -> pulsegradient_amplitude(t,j,pulse), t̄)
            G[j,q] = ∇Ω[:,q] · ∂Ω∂x
        end

        # SET THE GRADIENTS G[j] ~= ∂E/∂ν
        Δν = h * Utils.basisvector(length(x), (W+1)*q)
        G[W+1,q] = (f(x+Δν) - f0) / h
    end

end



##########################################################################################

# CONFIGURE BOUNDS
lbounds = Matrix{Float64}(undef, W+1, n); ubounds = Matrix{Float64}(undef, W+1, n)
# # #
lbounds[end,:] .= device.ω .- ΔΔ            # SET FREQUENCY CONSTRAINTS
ubounds[end,:] .= device.ω .+ ΔΔ
# # #
lbounds[1:W,:] .= -ΔΩ                       # SET MAXIMUM AMPLITUDE
ubounds[1:W,:] .=  ΔΩ

# CONFIGURE INITIAL PARAMETERS
Random.seed!(seed)
Ω0 = ΔΩ * (2*rand(W,n) .- 1)                # RANDOMLY SELECTED AMPLITUDES
ν0 = ΔΔ * (2*rand(n) .- 1) .+ device.ω      # RANDOMLY SELECTED FREQUENCIES
x0 = vec(vcat(Ω0, ν0'))

# OBTAIN INTIAL ENERGY ERROR
E0 = f(x0)
ε0 = E0 - E_FCI

# RUN L-BFGS-B OPTIMIZATION
_, x = lbfgsb(f, g!, x0, lb=vec(lbounds), ub=vec(ubounds),
    m=5,
    factr=1e4,      # UNITS OF MACHINE EPSILON: 1e9 ε ~= 1e-6
    pgtol=1e-5,
    iprint=50,
    maxfun=10000,
    maxiter=maxiter,
)

pulses = constructpulses(x)     # CONSTRUCT FINAL PULSES

##########################################################################################

# EXTRACT AMPLITUDES OVER THE PULSE DURATION
Ω = Matrix{Float64}(undef, r+1, n)      # CONSTRUCT Ω(t) FOR EACH QUBIT
for (q, pulse) in enumerate(pulses)
    Ω[:,q] .= map(t -> Pulses.amplitude(pulse, t), t̄)
end
ν = [Pulses.frequency(pulse, 0.0) for pulse in pulses]

# CALCULATE ENERGY AND OVERLAP OVER THE PULSE DURATION
E = Vector{Float64}(undef, length(t̄))
ε = Vector{Float64}(undef, length(t̄))
function callback(i, t, ψ)
    ψ = UD' * ψ                 # ROTATE INTO THE DEVICE BASIS
    ψ .*= exp.( (im*t) * ΛD)    # APPLY PHASES FOR INTERACTION PICTURE
    ψ = UD  * ψ                 # ROTATE OUT OF THE DEVICE BASIS
    E[i] = real(Utils.expectation(O, ψ))
    ε[i] = E[i] - E_FCI
end
ψ .= Evolutions.evolve(ψ0, pulses, device, Evolutions.Rotate;
        callback=callback, evolve_kwargs...)


# CALCULATE THE GRADIENT OVER THE PULSE DURATION
Gradients.gradientsignal(ψ0, pulses, device, O; gradient_kwargs...) # (assigns ∇Ω)

##########################################################################################

using FFTW
using FFTViews

# FOURIER TRANSFORM OUR TIME SERIES
Ω̂  = hcat((fft(Ω[:,q]) for q in 1:n)...)
∇Ω̂ = hcat((fft(∇Ω[:,q]) for q in 1:n)...)
Ê = fft(E)
ε̂ = fft(ε)

# NORMALIZE AND TRUNCATE OUT THE REDUNDANT NEGATIVE FREQUENCY COMPONENTS
ω̄ = 2π/T * (0:r÷2)              # ANGULAR FREQUENCIES
f̄ = ω̄ / 2π                      #  ACTUAL FREQUENCIES
Ω̂  =  Ω̂[1:(r÷2)+1,:] / (r+1)
∇Ω̂ = ∇Ω̂[1:(r÷2)+1,:] / (r+1)
Ê  =  Ê[1:(r÷2)+1]   / (r+1)
ε̂  =  ε̂[1:(r÷2)+1]   / (r+1)

##########################################################################################

# Ω(t) PLOT
plt_Ω = plot(
    ylabel="Amplitude (GHz)",
    ylim = [-ΔΩ, ΔΩ],
    legend=false,
    color_palette=:darkrainbow,
)
for q in 1:n
    plot!(plt_Ω, t̄, Ω[:,q]; color=q)
end

# ε(t) PLOT
plt_ε = plot(
    ylabel="Energy Error (Ha)",
    ylim = [0, max(ε...)],
    legend=false,
    color_palette=:darkrainbow,
)
plot!(plt_ε, t̄, ε; color=:black)

# ∇Ω(t) PLOT
plt_∇Ω = plot(
    xlabel="Time (ns)",
    ylabel="Gradient (Ha/GHz)",
    # ylim = [-ΔΩ, ΔΩ],
    legend=false,
    color_palette=:darkrainbow,
)
for q in 1:n
    plot!(plt_∇Ω, t̄, ∇Ω[:,q]; color=q)
end




# INDEX TO TRUNCATE FREQUENCY PLOT
maxω = 44
k = findfirst(ω -> ω > maxω, ω̄)

# Ω(ω) PLOT
plt_Ω̂ = plot(
    ylabel="Amplitude (GHz)",
    legend=false,
    color_palette=:darkrainbow,
)
for q in 1:n
    plot!(plt_Ω̂, ω̄[2:k], abs.(Ω̂[2:k,q]); color=q)
end

# ε(ω) PLOT
plt_ε̂ = plot(
    ylabel="Energy Error (Ha)",
    legend=false,
    color_palette=:darkrainbow,
)
plot!(plt_ε̂, ω̄[2:k], abs.(ε̂[2:k,1]); color=:black)

# ∇Ω(ω) PLOT
plt_∇Ω̂ = plot(
    xlabel="Frequency (GHz)",
    ylabel="Gradient (Ha/GHz)",
    legend=false,
    color_palette=:darkrainbow,
)
for q in 1:n
    plot!(plt_∇Ω̂, ω̄[2:k], abs.(∇Ω̂[2:k,q]); color=q)
end


# ∇Ω(ω) PLOT WITH DETUNING BARS
plt_∇Ω̂_ν = plot(
    xlabel="Frequency (GHz)",
    ylabel="Gradient (Ha/GHz)",
    legend=false,
    color_palette=:darkrainbow,
)
for q in 1:n
    plot!(plt_∇Ω̂_ν, ω̄[2:k], abs.(∇Ω̂[2:k,q]); color=q)
end

for q in 1:n
    vline!(plt_∇Ω̂_ν, [abs(ν0[q]-device.ω[q])], color=q, linestyle=:dot)
    vline!(plt_∇Ω̂_ν, [abs( ν[q]-device.ω[q])], color=q, linestyle=:dash)
end


##########################################################################################

function tickformatter(x)
    #= Objective: return a string with *exactly* 10 characters. =#
    maxlength = 10
    ecutoff = 5

    if x == 0
        return " "^(maxlength-1) * "0"
    end

    power = findlast(p -> 10.0^p <= abs(x), (-15):15) - 16
    mantissa = abs(x) / 10.0^power


    if abs(power) >= ecutoff
        # USE SCIENTIFIC NOTATION
        signstring = sign(x) == -1 ? "-" : " "
        powerstring = "e$power"
        digitsleft = maxlength - 1 - length(powerstring)
        mantissastring = string(round(mantissa, digits=digitsleft-2))
        numstring = signstring * mantissastring * powerstring
    else
        # WRITE THE NUMBER AS IT IS
        numstring = string(x)
        if length(numstring) > maxlength
            numstring = SubString(numstring, 1, maxlength)
        end

        # REMOVE TRAILING ZEROS
        if '.' in numstring
            while numstring[end] == '0'
                numstring = numstring[begin:end-1]
            end
        end
    end

    if length(numstring) < maxlength
        numstring = " "^(maxlength-length(numstring)) * numstring
    end

    return numstring

end

# ALL TIME SERIES
plot(plt_Ω, plt_ε, plt_∇Ω,
    legend=false,
    layout=(3,1),
    dpi=300,
    yformatter=tickformatter,
    tickfont=font("courier bold", 4),
    guidefont=font("courier bold", 6),
)
savefig("fig/cambridge.timeseries.$maxiter.png")
savefig("fig/cambridge.timeseries.$maxiter.pdf")

# ALL SPECTRA
plot(plt_Ω̂, plt_ε̂, plt_∇Ω̂,
    legend=false,
    layout=(3,1),
    dpi=300,
    guidefontsize=8,
    yformatter=tickformatter,
    tickfont=font("courier bold", 4),
    guidefont=font("courier bold", 6),
)
savefig("fig/cambridge.spectra.$maxiter.png")
savefig("fig/cambridge.spectra.$maxiter.pdf")

# ALL SPECTRA, WITH DETUNING BARS
plot(plt_Ω̂, plt_ε̂, plt_∇Ω̂_ν,
    legend=false,
    layout=(3,1),
    dpi=300,
    guidefontsize=8,
    yformatter=tickformatter,
    tickfont=font("courier bold", 4),
    guidefont=font("courier bold", 6),
)
savefig("fig/cambridge.spectrabars.$maxiter.png")
savefig("fig/cambridge.spectrabars.$maxiter.pdf")

