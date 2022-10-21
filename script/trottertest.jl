#=
As a brief Julia tutorial,
    let us quantify the Trotter error associated with Ω=exp[𝑖Δ (z â + h.t.)]

The idea is to decompose z = x + 𝑖y, and Ω as exp[𝑖Δ (x Q̂ + y P̂)],
    then approximate this as exp(𝑖Δ x Q̂) exp(𝑖Δ y P̂),
    each of which can be pre-diagonalized.

=#

using LinearAlgebra
using Random

# DEFINE THE NUMBER OF STATES
n = 3

# CONSTRUCT THE ANNIHILATION OPERATOR (HARMONIC OSCILLATOR NORMALIZATION)
a = zeros((n,n))
for i ∈ 1:n-1
    a[i,i+1] = √i
end

# CONSTRUCT CANONICAL OPERATORS
Q = a + a'
P = im*(a - a')

# DIAGONALIZE CANONICAL OPERATORS
ΛQ, UQ = eigen(Q)
ΛP, UP = eigen(P)
UQP = UQ' * UP
UPQ = UP' * UQ

# DEFINE FUNCTION TO PERFORM EXPONENTIATION
function exponentiate(x, y, Δ=1.0, p=0)
    if p == 0                   # EXACT EXPONENTIATION
        T = complex(x,y) * a    # CLUSTER OPERATOR, OF SORTS
        exp(im * Δ * (T + T'))
    elseif p == 1               # PRODUCT FORMULA
        exp(im*Δ*x*Q) * exp(im*Δ*y*P)
    elseif p == 2               # SYMMETRIZED PRODUCT FORMULA
        exp(im*Δ*x/2*Q) * exp(im*Δ*y*P) * exp(im*Δ*x/2*Q)
    end
end

# DEFINE FUNCTION TO PERFORM EXPONENTIATION
function prediagonalizedexponentiate(x, y, Δ=1.0, p=0)
    if p == 0                   # EXACT EXPONENTIATION
        T = complex(x,y) * a    # CLUSTER OPERATOR, OF SORTS
        exp(im * Δ * (T + T'))
    elseif p == 1               # PRODUCT FORMULA
        expQ = exp.(im*Δ*x*ΛQ)
        expP = exp.(im*Δ*y*ΛQ)
        UQ * Diagonal(expQ) * UQP * Diagonal(expP) * UP'
    elseif p == 2               # SYMMETRIZED PRODUCT FORMULA
        expQ = exp.(im*Δ*x/2*ΛQ)    # HALVED SO THAT IT FLANKS P FACTOR
        expP = exp.(im*Δ*y*ΛQ)
        UQ * Diagonal(expQ) * UQP * Diagonal(expP) * UPQ * Diagonal(expQ) * UQ'
    end
end

# DEFINE FUNCTION TO MEASURE DIFFERENCE BETWEEN SOLUTIONS
function distance(A, B)
    # I'm not really sure what to use here...
    opnorm(A-B)
end

# DEFINE EXPERIMENT
function experiment(x, y, Δ)
    exact = exponentiate(x,y,Δ)
    ε1 = distance(exponentiate(x,y,Δ,1), exact)
    ε2 = distance(exponentiate(x,y,Δ,2), exact)
    ε1pd = distance(prediagonalizedexponentiate(x,y,Δ,1), exact)
    ε2pd = distance(prediagonalizedexponentiate(x,y,Δ,2), exact)
    return ε1, ε2, ε1pd, ε2pd
end

# GENERATE STATISTICS
N = 1000                    # NUMBER OF TRIALS
x = rand(Float64, N)
y = rand(Float64, N)
logΔ = -6*rand(Float64, N)

ε = Array{Float64}(undef, (N,4))
for i in 1:N
    ε[i,:] .= experiment(x[i], y[i], 10^logΔ[i])
end

# DISPLAY RESULTS
# import Pkg
# Pkg.add("Plots")
using Plots
plt = scatter(logΔ, log10.(ε[:,1]))
scatter!(logΔ, log10.(ε[:,2]))
scatter!(logΔ, log10.(ε[:,3]))
scatter!(logΔ, log10.(ε[:,4]))

# gui(plt)
