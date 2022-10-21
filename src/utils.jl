#= Provide matrix representations for common Hamiltonians. =#
module Utils
export a_matrix, on

using LinearAlgebra


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

end
