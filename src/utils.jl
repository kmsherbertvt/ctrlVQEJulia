#= Provide matrix representations for common Hamiltonians. =#
module Utils
export a_matrix, on

import LinearAlgebra: kron, I


"""
    a_matrix(nstates::Integer=2)

Matrix representation of the annihilation operator on a single qubit.

This method assumes a harmonic system truncated to `nstates` levels.

TODO: Is this generally useful for qudits, or just transmons?

"""
function a_matrix(nstates::Integer=2)
    a = zeros((nstates,nstates))
    for i ∈ 1:nstates-1
        a[i,i+1] = √i               # √i ENSURES STATES PRESERVE NORMALIZATION
    end                             # TODO: Kyle, please re-do the math to understand why.
    return a
end

"""
    on(op::AbstractMatrix{<:Number}, q::Integer, n::Integer)

Expand the single-qubit matrix operator `op` to act on qubit `q` of `n` qubits.

In other words, apply Kronecker-products such that identity `I` acts on the other qubits.

"""
function on(op::AbstractMatrix{<:Number}, q::Integer, n::Integer)
    A = ones(1,1)                   # A 1x1 IDENTITY MATRIX
    I = one(op)                     # AN IDENTITY MATRIX MATCHING THE DIMENSIONS OF `op`
    for i ∈ 1:n
        A = kron(A, i == q ? op : I)    # `op` ACTS ON QUBIT `q`, `I` ACTS ON ALL OTHERS
    end
    return A
end

"""
    kron_concat(ops::AbstractVector{<:AbstractMatrix{<:Number}})

Concatenate a sequence of operators with the Kronecker product.

"""
function kron_concat(ops::AbstractVector{<:AbstractMatrix{<:Number}})
    O = Matrix(I,1,1)
    for q ∈ eachindex(ops)
        O = kron(O, ops[q])
    end
    return O
    #= TODO: Pre-allocate O, ...and each intermediate. :?
        There *must* be an easy way to kron a sequence of matrices together...
        Maybe reshaping O into a tensor..? Ironic...
    =#
end

"""
    kron_concat(ops::AbstractMatrix{<:Number}, n::Integer)

Concatenate a repeated string of an operator with the Kronecker product.

"""
function kron_concat(op::AbstractMatrix{<:Number}, n::Integer)
    O = Matrix(I,1,1)
    for q ∈ 1:n
        O = kron(O, op)
    end
    return O
    #= TODO: Pre-allocate O, ...and each intermediate. :?
        There *must* be an easy way to kron a sequence of matrices together...
        Maybe reshaping O into a tensor..? Ironic...
    =#
end



"""
    algebra(
        n::Integer,
        m::Integer=2;
        basis::Union{AbstractMatrix{<:Number},Nothing}=nothing,
    )

Construct a vector of annihilation operators acting on each of n m-level systems.

Optionally, rotate these operators into the provided basis.

These matrices, in conjunction with their adjoints,
    form a complete algebra for the space of n m-level systems.

"""
function algebra(
    n::Integer,
    m::Integer=2;
    basis::Union{AbstractMatrix{<:Number},Nothing}=nothing,
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


end
