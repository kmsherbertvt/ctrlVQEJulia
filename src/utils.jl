#= Provide matrix representations for common Hamiltonians. =#

module Utils

import LinearAlgebra: kron, I, eigen, Eigen, Hermitian, norm, mul!

"""
    basisvector(N, i)

Construct the standard basis vector e⃗ᵢ.

"""
function basisvector(N, i)
    e = zeros(Bool, N)
    e[i] = 1
    return e
end

"""
    transform!(x::Vector{T}, A::Matrix{<:Number}[, _x::Vector{T}]) where T <: Number

Calculate the matrix-vector product A·x and copy the result into `x`.

`_x` is a pre-allocated "work variable" storing intermediate results.
In fact, this method just calls Julia's in-place matrix multiplication,
    then copies the results to the original vector.
"""
function transform!(
    x::Vector{T}, A::AbstractMatrix{<:Number}, _x::Vector{T}
) where T <: Number
    x .= mul!(_x, A, x)
end

transform!(x::Vector{<:Number}, A::AbstractMatrix{<:Number}) = transform!(
    x, A, Vector{eltype(x)}(undef, length(x))
)

"""
    braket(
        φ::Vector{<:Number}, A::Matrix{<:Number}, ψ::Vector{T}[, _ψ::Vector{T}]
    ) where T <: Number

Calculate the scalar quantity ⟨φ|A|ψ⟩.

`_ψ` is a pre-allocated "work variable" storing the intermediate result A|ψ⟩,
    so its dimension should match that of ψ.

"""
function braket(
    φ::Vector{<:Number}, A::Matrix{<:Number}, ψ::Vector{T}, _ψ::Vector{T}
) where T <: Number
    mul!(_ψ, A, ψ)
    return σ' * _ψ
end

braket(φ::Vector{<:Number}, A::Matrix{<:Number}, ψ::Vector{<:Number}) = braket(
    φ, A, ψ, Vector{eltype(ψ)}(undef, length(ψ))
)

"""
    on(op::Matrix{<:Number}, q::Int, n::Int)

Expand the single-qubit matrix operator `op` to act on qubit `q` of `n` qubits.

In other words, apply Kronecker-products such that identity `I` acts on the other qubits.

"""
function on(op::Matrix{<:Number}, q::Integer, n::Integer)
    A = ones(1,1)                   # A 1x1 IDENTITY MATRIX
    I = one(op)                     # AN IDENTITY MATRIX MATCHING THE DIMENSIONS OF `op`
    for i ∈ 1:n
        A = kron(A, i == q ? op : I)    # `op` ACTS ON QUBIT `q`, `I` ACTS ON ALL OTHERS
    end
    return A
end

"""
    kron_concat!(
        ops::AbstractVector{Matrix{T}},
        [O_::AbstractVector{Matrix{T}},]
    ) where T <: Number

Concatenate a sequence of operators with the Kronecker product.

If `O_` is provided, each successive pairwise kron operation
    is written to a pre-allocated element of O_.

"""
function kron_concat(
    ops::AbstractVector{Matrix{T}},
    O_::AbstractVector{Matrix{T}},
) where T <: Number
    O_[1] .= ops[1]
    for q ∈ 2:length(ops)
        kron!(O_[q], O_[q-1], ops[q])
    end
    return O_[end]
end

function kron_concat(ops::AbstractVector{<:Matrix{<:Number}})
    O = Matrix(I,1,1)
    for q ∈ eachindex(ops)
        O = kron(O, ops[q])
    end
    return O
end

"""
    kron_concat(
        op::Matrix{T}, n::Integer,
        [O_::AbstractVector{Matrix{T}},]
    ) where T <: Number

Concatenate a repeated string of an operator with the Kronecker product.

If `O_` is provided, each successive pairwise kron operation
    is written to a pre-allocated element of O_.

"""
function kron_concat(
    op::Matrix{T}, n::Integer,
    O_::AbstractVector{Matrix{T}}
) where T <: Number
    O_[1] .= op
    for q ∈ 2:n
        kron!(O_[q], O_[q-1], op)
    end
    return O_[end]
end

function kron_concat(op::Matrix{<:Number}, n::Integer)
    O = Matrix(I,1,1)
    for q ∈ 1:n
        O = kron(O, op)
    end
    return O
end


"""
    qubitspace(
        ψ::Vector{ComplexF64},              # `m`-LEVEL STATEVECTOR
        n::Integer;                         # NUMBER OF QUBITS

        # INFERRED VALUES (relatively fast, but pass them in to minimize allocations)
        N = length(ψ),                      # SIZE OF STATEVECTOR
        m = round(Int, N^(1/n)),            # NUMBER OF LEVELS ON EACH QUBIT

        # PRE-ALLOCATIONS (for those that want every last drop of efficiency...)
        ψ2= Vector{ComplexF64}(undef, 2^n), # 2-LEVEL STATEVECTOR
        z = Vector{Bool}(undef, n),         # BITSTRING VECTOR
    )

Project an `m`-level system of `n` qubits onto a 2-level system of `n` qubits.

The resulting statevector is *not* normalized (you'll probably want to do so).

The value of `m` is inferred from the length of `ψ`.

"""
function qubitspace(
    ψ::Vector{T},                       # `m`-LEVEL STATEVECTOR
    n::Integer;                         # NUMBER OF QUBITS

    # INFERRED VALUES (relatively fast, but pass them in to minimize allocations)
    N = length(ψ),                      # SIZE OF STATEVECTOR
    m = round(Int, N^(1/n)),            # NUMBER OF LEVELS ON EACH QUBIT

    # PRE-ALLOCATIONS (for those that want every last drop of efficiency...)
    ψ2= Vector{T}(undef, 2^n),          # 2-LEVEL STATEVECTOR
    z = Vector{Bool}(undef, n),         # BITSTRING VECTOR
) where T <: Number
    # SELECT ELEMENTS OF ψ WITH ONLY 0, 1 VALUES
    for i2 in eachindex(ψ2)
        digits!(z, i2-1, base=2)                        # FILL z WITH i2's BITSTRING
        i = 1+foldr((a,b) -> muladd(m,b,a), z, init=0)  # PARSE z AS BASE `m` STRING
        ψ2[i2] = ψ[i]
    end

    return ψ2
end

#= TODO: We can generally transmute any m-level statevector to an m'-level statevector.

`qubitspace` method is easily generalized to arbitrary m' when m > m'. Unnormalized.
When m' > m, start with all zeros, and then copy in values you have. Normalized.
    Latter is basically a 1D version of `extendoperator` below.

=#

function extendoperator(
    A::AbstractMatrix,      # OPERATOR
    n::Integer,             # NUMBER OF QUBITS
    m::Integer;             # TARGET NUMBER OF LEVELS ON EACH QUBIT

    # INFERRED VALUES (relatively fast, but pass them in to minimize allocations)
    N = m^n,                            # SIZE OF EXTENDED HILBERT SPACE
    N0 = size(A,1),                     # SIZE OF INITIAL HILBERT SPACE
    m0 = round(Int, N0^(1/n)),          # INITIAL NUMBER OF LEVELS ON EACH QUBIT

    # PRE-ALLOCATIONS (for those that want every last drop of efficiency...)
    B = zeros(eltype(A), (N,N)),        # EXTENDED OPERATOR (STARTS WITH ZEROS!)
    imap = nothing,                     # MAP FROM BASE-m0 INDICES TO BASE-m INDICES
)
    if imap === nothing
        z = Vector{Int}(undef, n)       # BITSTRING VECTOR
        imap = Vector{Int}(undef, N0)   # MAP FROM BASE-m0 INDICES TO BASE-m INDICES
        for i in 1:N0
            digits!(z, i-1, base=m0)                    # FILL z WITH i's BASE `m0` STRING
            imap[i] = 1+foldr((a,b)->muladd(m,b,a), z, init=0)# PARSE z AS BASE `m` STRING
        end
    end

    # COPY VALUES OF A INTO B
    for i in 1:N0; for j in 1:N0
        B[imap[i],imap[j]] = A[i,j]
    end; end

    return B
end

#= TODO: I think there should be a well-defined `shrinkoperator` method also. =#


"""
    projector(n::Integer, m::Integer, m0::Integer)

Project a Hilbert space of `n` `m0`-level qubits onto that of `n` `m`-level qubits

Returns an (`n^m`, `n^m0`) shaped matrix `Π`.
To perform the projection on a vector, use ψ ← Πψ.
To perform the projection on a matrix, use A ← ΠAΠ'.

"""
function projector(n::Integer, m::Integer, m0::Integer)
    # TODO: Pre-allocations
    # IT'S EASIER TO CALCULATE PROJECTOR FROM SMALLER SPACE TO LARGER SPACE
    if m < m0; return projector(n, m0, m)'; end     # NOTE THE ADJOINT ' OPERATOR

    z = Vector{Int}(undef, n)       # PRE-ALLOCATION TO STORE BASE-m0 DECOMPOSITIONS
    N  = m^n; N0 = m0^n             # FULL HILBERT SPACE DIMENSIONS
    Id= Matrix{Bool}(I, N, N)       # IDENTITY MATRIX IN LARGER HILBERT SPACE
    Π = Matrix{Bool}(undef, N, N0)  # PROJECTOR MATRIX
    j = 1                           # ITERATES COLUMNS OF PROJECTOR

    for i ∈ 1:N                     # FILL PROJECTOR WITH SELECT COLUMNS FROM IDENTITY
        digits!(z, i-1, base=m)         # DECOMPOSE INDEX INTO LARGER BASE
        if any(z .>= m0); continue; end # SKIP INDEX IF ABSENT IN SMALLER SPACE
        Π[:,j] .= Id[:,i]               # COPY COLUMN OF IDENTITY MATRIX
        j += 1                          # ADVANCE COLUMN INDEX
    end

    # for i ∈ 1:N                     # FILL PROJECTOR WITH SELECT COLUMNS FROM IDENTITY
    #     digits!(z, i-1, base=m)         # DECOMPOSE INDEX INTO LARGER BASE
    #     if any(z .>= m0); continue; end # SKIP INDEX IF ABSENT IN SMALLER SPACE
    #     j = 1+foldr((a,b)->muladd(m0,b,a), reverse(z), init=0)
    #     Π[:,j] .= Id[:,i]               # COPY COLUMN OF IDENTITY MATRIX
    # end


    return Π
end






"""
    expectation(A::Matrix{ComplexF64}, ψ::Vector{ComplexF64})

Evaluate the expectation value ⟨ψ|A|ψ⟩.

`A` and `ψ` should have compatible dimensions.

"""
function expectation(
    A::AbstractMatrix{<:Number},
    ψ::Vector{T};

    # INFERRED VALUES (relatively fast, but pass them in to minimize allocations)
    N = length(ψ),                          # SIZE OF STATEVECTOR

    # PRE-ALLOCATIONS (for those that want every last drop of efficiency...)
    tmpV = Vector{T}(undef, N),    # FOR MATRIX-VECTOR MULTIPLICATION
) where T <: Number
    mul!(tmpV, A, ψ)
    return ψ' * tmpV
end

"""
    a_matrix(m::Int=2)

Matrix representation of the bosonic annihilation operator on a single qubit.

Note that the matrix representation must be truncated to `m` levels.

"""
function a_matrix(m::Integer=2)
    a = zeros((m,m))
    for i ∈ 1:m-1
        a[i,i+1] = √i               # BOSONIC ANNIHILATION OPERATOR
    end
    return a
end

"""
    algebra(
    n::Int,
    m::Int=2;
    basis::Union{Matrix{<:Number},Nothing}=nothing,
)

Construct a vector of annihilation operators acting on each of `n` `m`-level systems.

Optionally, rotate these operators into the provided basis.

These matrices, in conjunction with their adjoints,
    form a complete algebra for the space of `n` `m`-level systems.

"""
function algebra(
    n::Int,
    m::Int=2;
    basis::Union{Matrix{<:Number},Nothing}=nothing,
)
    a_ = a_matrix(m)                        # SINGLE-QUBIT ANNIHILATION OPERATOR
    a = [on(a_, q, n) for q in 1:n]         # EACH OPERATOR, ACTING ON FULL HILBERT SPACE
    if !(basis === nothing)
        for q ∈ 1:n
            a[q] .= basis' * a[q] * basis   # CONJUGATE WITH BASIS
        end
    end
    return a
end


"""
    dressedbasis(H::Hermitian)

Calculate the eigenvalues and eigenvectors for a Hermitian matrix `H`.

This method returns an `Eigen` object, the same type as LinearAlgebra.eigen(H).

Dressing does three things:
1. Re-order eigenpairs such that the ith eigenvector is the one with
    the largest magnitude in the ith element.
    (This is distinct from the LinearAlgebra.eigen order,
        which uses ascending order of eigenvalues.)
2. Multiply each eigenvector by a global phase such that the ith element
    is real and positive.
3. Zero out any elements with magnitude less than machine precision.

"""
function dressedbasis(H::Hermitian)
    N = size(H)[1]
    Λ, U = eigen(H)

    # IMPOSE PERMUTATION
    σ = Vector{Int}(undef,N)
    for i in 1:N
        perm = sortperm(abs.(U[i,:]), rev=true) # STABLE SORT BY "ith" COMPONENT
        perm_ix = 1                             # CAREFULLY HANDLE TIES
        while perm[perm_ix] ∈ σ[1:i-1]
            perm_ix += 1
        end
        σ[i] = perm[perm_ix]
    end
    Λ .= Λ[  σ]
    U .= U[:,σ]

    #= TODO: This strategy for selecting the permutation might not be the most efficient.
        It's not *bad*: should be no worse than O(N² log N).
        But a less careful handling of ties would need only O(N²).
        So I feel like there ought be a strategy that doesn't involve a full sort...
    =#

    # IMPOSE PHASE
    for i in 1:N
        U[:,i] .*= U[i,i] < 0 ? -1 : 1           # ASSUMES REAL TYPE
        # U[:,i] .*= exp(-im*angle(U[i,i]))        # ASSUMES COMPLEX TYPE
    end

    # IMPOSE ZEROS
    Λ[abs.(Λ) .< eps(eltype(Λ))] .= zero(eltype(Λ))
    U[abs.(U) .< eps(eltype(U))] .= zero(eltype(U))
    return Eigen(Λ,U)

    #= TODO: Each of these impositions should properly be their own function.
        The phase and zeros impositions should have two methods each, Real vs Complex types.

            real phase should just multiply U[:,i] by -1 if U[i,i] < 0.
            complex phase does as written
                except use .*= operator once you're guaranteed types needn't change.

            complex zeros should reinterpret as (larger) Float array
                and zero out all small components.
            I don't know exactly how to do the reinterpret with the correct typing,
                short of writing a different method for each possible complex type.
            I absolutely do *not* want to assume Float64.
                It's *probably* true but I don't trust it!

            phase need only ever take matrices, but zeros takes anything
            Unfortunately, knowing the NUMBER type probably means
                having to dispatch separately on the array type.
            I think we can do that just fine with a Union rather than separate methods,
                but you can see how this whole thing gets WAY more complicated
                than it should be when you want to zero out individual components.
            abs turns out to be WAY easier...
    =#
end

"""
    infidelity(ψ0, ψ)

Calculates ``1 - |⟨ψ0|ψ⟩|^2``.

"""
infidelity(ψ,φ) = 1 - abs2(ψ'*φ)



"""
    heaviside(t)

The heaviside step function:

            { 0     t < 0
    Θ(t) =  { 1/2   t = 0
            { 1     t > 0

"""
function heaviside(t)
    return (sign(t) + 1) / 2
end

"""
    interval(t, a, b)

An interval function, which gives 0 when t is outside [a,b] and 1 inside.

This function is implemented with heaviside step-functions,
    so that 1/2 is returned when t is on the bounds.

"""
function interval(t, a, b)
    return heaviside(t-a) - heaviside(t-b)
end

"""
    diracdelta(t, τ)

An approximation of the Dirac delta distribution, under a finite step-size τ.

The essential property of the Dirac delta is:

    ∫dt δ(t) = 1

    when the bounds of the integral include 0.

Under a finite step-size, this changes to:

    ∑ τ·δ(t) = 1

    when the bounds of the sum would include t=0.

"""
function diracdelta(t, τ)
    return interval(t, -τ/2, τ/2) / τ
end


end # END MODULE
